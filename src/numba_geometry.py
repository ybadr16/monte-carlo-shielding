"""Flat-array CSG geometry for the experimental Numba engine.

The object-based geometry (Plane/Sphere/Cylinder/Box/Region) is compiled once
into flat NumPy arrays that ``@njit`` kernels can traverse without Python
objects: a surface table (type code + coefficients + boundary condition) and,
per medium, a postfix (RPN) token stream encoding its Boolean CSG expression.
The njit stack evaluator walks the postfix to answer ``contains`` for arbitrary
nested intersection/union/complement/difference regions, matching
medium.Region.contains, and ``nearest_crossing`` mirrors
vector_transport.nearest_crossings (nearest primitive whose crossing point lies
in its medium). This is the load-bearing piece for a fully-njit transport
engine; everything downstream is data-flattening on top of it.
"""
import math

import numpy as np

from .medium import Plane, Sphere, Cylinder, Region

try:
    from numba import njit
    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False

    def njit(*a, **k):
        return (lambda f: f) if not a else a[0]

# surface type codes
PLANE, SPHERE, CYLINDER = 0, 1, 2
# postfix operator tokens (surface tests are encoded as token = surface_index)
OP_AND, OP_OR, OP_NOT = -1, -2, -3
_TOL = 1e-9


# ---------------------------------------------------------------------------
# njit surface primitives (mirror medium.py's evaluate / nearest_surface_method)
# ---------------------------------------------------------------------------

@njit(fastmath=False, cache=True)
def surf_evaluate(stype, p, x, y, z):
    """Signed surface function; <= tolerance means the point is 'inside'."""
    if stype == PLANE:
        return p[0] * x + p[1] * y + p[2] * z + p[3]
    if stype == SPHERE:
        dx = x - p[0]; dy = y - p[1]; dz = z - p[2]
        return dx * dx + dy * dy + dz * dz - p[3] * p[3]
    # cylinder: p[4] axis 0=x,1=y,2=z
    ax = p[4]
    if ax == 0.0:
        a = y - p[1]; b = z - p[2]
    elif ax == 1.0:
        a = x - p[0]; b = z - p[2]
    else:
        a = x - p[0]; b = y - p[1]
    return a * a + b * b - p[3] * p[3]


@njit(fastmath=False, cache=True)
def surf_distance(stype, p, x, y, z, u, v, w):
    """Forward distance to the surface along (u,v,w); np.inf if none.
    Branch-for-branch mirror of the scalar nearest_surface_method."""
    if stype == PLANE:
        num = -p[3] - p[0] * x - p[1] * y - p[2] * z
        den = p[0] * u + p[1] * v + p[2] * w
        if abs(num) <= 1e-8:
            return 0.0 if abs(den) > 1e-8 else np.inf
        if den == 0.0:
            return np.inf
        t = num / den
        return t if t >= 0.0 else np.inf
    if stype == SPHERE:
        dn = math.sqrt(u * u + v * v + w * w)
        if dn == 0.0:
            return np.inf
        un = u / dn; vn = v / dn; wn = w / dn
        xb = x - p[0]; yb = y - p[1]; zb = z - p[2]
        k = xb * un + yb * vn + zb * wn
        c = xb * xb + yb * yb + zb * zb - p[3] * p[3]
        disc = k * k - c
        if disc < 0.0:
            return np.inf
        sq = math.sqrt(disc); d1 = -k - sq; d2 = -k + sq
        if c < 0.0:
            return d1 if d1 > d2 else d2
        if d1 >= 0.0:
            return d1
        return d2 if d2 >= 0.0 else np.inf
    # cylinder
    ax = p[4]
    if ax == 0.0:      # x-axis: use (y,z) / (v,w)
        rb1 = y - p[1]; rb2 = z - p[2]; a = v * v + w * w; k = rb1 * v + rb2 * w; d1c = v; d2c = w
    elif ax == 1.0:    # y-axis: use (x,z) / (u,w)
        rb1 = x - p[0]; rb2 = z - p[2]; a = u * u + w * w; k = rb1 * u + rb2 * w; d1c = u; d2c = w
    else:              # z-axis: use (x,y) / (u,v)
        rb1 = x - p[0]; rb2 = y - p[1]; a = u * u + v * v; k = rb1 * u + rb2 * v; d1c = u; d2c = v
    c = rb1 * rb1 + rb2 * rb2 - p[3] * p[3]
    if a == 0.0:
        return np.inf
    disc = k * k - a * c
    if disc < 0.0:
        return np.inf
    sq = math.sqrt(disc); d1 = (-k - sq) / a; d2 = (-k + sq) / a
    if c < 0.0:
        return d1 if d1 > d2 else d2
    if c == 0.0:
        dot = rb1 * d1c + rb2 * d2c
        if dot < 0.0 and d2 >= 0.0:
            return d2
        return np.inf
    lo = np.inf
    if d1 >= 0.0 and d1 < lo:
        lo = d1
    if d2 >= 0.0 and d2 < lo:
        lo = d2
    return lo


@njit(cache=True)
def contains_medium(m, tok, tok_off, stype, params, x, y, z):
    """Evaluate medium ``m``'s postfix CSG expression at a point.
    Boolean stack machine over the token stream; matches Region.contains."""
    stack = np.empty(32, dtype=np.uint8)
    sp = 0
    for t in range(tok_off[m], tok_off[m + 1]):
        tk = tok[t]
        if tk >= 0:
            inside = surf_evaluate(stype[tk], params[tk], x, y, z) <= _TOL
            stack[sp] = 1 if inside else 0
            sp += 1
        elif tk == OP_AND:
            b = stack[sp - 1]; a = stack[sp - 2]; sp -= 1
            stack[sp - 1] = 1 if (a and b) else 0
        elif tk == OP_OR:
            b = stack[sp - 1]; a = stack[sp - 2]; sp -= 1
            stack[sp - 1] = 1 if (a or b) else 0
        else:  # OP_NOT
            stack[sp - 1] = 0 if stack[sp - 1] else 1
    return stack[0] == 1


@njit(cache=True)
def resolve_medium(n_media, priority, tok, tok_off, stype, params, x, y, z):
    """Highest-priority medium containing the point, or -1 (strict >, so an
    earlier medium wins a priority tie, matching the scalar loop)."""
    best = -1
    best_pri = -1e30
    for m in range(n_media):
        if priority[m] > best_pri and contains_medium(m, tok, tok_off, stype, params, x, y, z):
            best = m; best_pri = priority[m]
    return best


@njit(cache=True)
def nearest_crossing(n_media, prim, prim_off, priority, tok, tok_off,
                     stype, params, x, y, z, u, v, w):
    """Nearest boundary crossing over all media: the smallest forward distance
    to a primitive whose crossing point lies inside that primitive's medium.
    Returns (distance, surface_index); (inf, -1) on escape. Mirrors
    vector_transport.nearest_crossings including region/surface iteration order.
    """
    best = np.inf
    best_s = -1
    for m in range(n_media):
        for pi in range(prim_off[m], prim_off[m + 1]):
            s = prim[pi]
            d = surf_distance(stype[s], params[s], x, y, z, u, v, w)
            if d >= 0.0 and d < best:
                px = x + d * u; py = y + d * v; pz = z + d * w
                if contains_medium(m, tok, tok_off, stype, params, px, py, pz):
                    best = d; best_s = s
    return best, best_s


# ---------------------------------------------------------------------------
# Python compiler: object geometry -> flat arrays
# ---------------------------------------------------------------------------

def _surf_row(surf):
    """(type_code, 6-param row, bc_code) for a primitive surface."""
    bc = 1 if getattr(surf, "boundary_type", "transmission") == "reflective" else 0
    if isinstance(surf, Plane):
        return PLANE, [surf.A, surf.B, surf.C, surf.D, 0.0, 0.0], bc
    if isinstance(surf, Sphere):
        return SPHERE, [surf.x0, surf.y0, surf.z0, surf.radius, 0.0, 0.0], bc
    if isinstance(surf, Cylinder):
        ax = {"x": 0.0, "y": 1.0, "z": 2.0}[surf.axis]
        return CYLINDER, [surf.x0, surf.y0, surf.z0, surf.radius, ax, 0.0], bc
    raise TypeError(f"unsupported surface {type(surf).__name__}")


class CompiledGeometry:
    """Flat-array form of a list of media, ready for the njit kernels."""

    def __init__(self, mediums):
        stype, params, bc = [], [], []
        surf_id = {}  # id(surface) -> flat index (dedupe shared surfaces)

        def intern(surf):
            key = id(surf)
            if key not in surf_id:
                t, row, b = _surf_row(surf)
                surf_id[key] = len(stype)
                stype.append(t); params.append(row); bc.append(b)
            return surf_id[key]

        def emit_postfix(region, out):
            """Append region's Boolean expression in postfix to ``out``."""
            op = region.operation
            surfs = region.surfaces

            def operand(s):
                if isinstance(s, Region):
                    emit_postfix(s, out)
                else:
                    out.append(intern(s))

            if op == "intersection":
                operand(surfs[0])
                for s in surfs[1:]:
                    operand(s); out.append(OP_AND)
            elif op == "union":
                operand(surfs[0])
                for s in surfs[1:]:
                    operand(s); out.append(OP_OR)
            elif op == "complement":
                operand(surfs[0]); out.append(OP_NOT)
            elif op == "difference":
                a, b = surfs
                operand(a); operand(b); out.append(OP_NOT); out.append(OP_AND)
            else:
                raise ValueError(f"unknown op {op}")

        def prims_of(region, out):
            for s in region.surfaces:
                if isinstance(s, Region):
                    prims_of(s, out)
                else:
                    out.append(intern(s))

        tok, tok_off = [], [0]
        prim, prim_off = [], [0]
        priority, is_void = [], []
        for med in mediums:
            emit_postfix(med, tok); tok_off.append(len(tok))
            prims_of(med, prim); prim_off.append(len(prim))
            priority.append(float(med.priority))
            is_void.append(1 if med.is_void else 0)

        self.n_media = len(mediums)
        self.stype = np.array(stype, dtype=np.int64)
        self.params = np.array(params, dtype=np.float64) if params else np.zeros((0, 6))
        self.bc = np.array(bc, dtype=np.int64)
        self.tok = np.array(tok, dtype=np.int64)
        self.tok_off = np.array(tok_off, dtype=np.int64)
        self.prim = np.array(prim, dtype=np.int64)
        self.prim_off = np.array(prim_off, dtype=np.int64)
        self.priority = np.array(priority, dtype=np.float64)
        self.is_void = np.array(is_void, dtype=np.int64)

    def contains(self, m, x, y, z):
        return contains_medium(m, self.tok, self.tok_off, self.stype, self.params, x, y, z)

    def resolve(self, x, y, z):
        return resolve_medium(self.n_media, self.priority, self.tok, self.tok_off,
                              self.stype, self.params, x, y, z)

    def crossing(self, x, y, z, u, v, w):
        return nearest_crossing(self.n_media, self.prim, self.prim_off, self.priority,
                                self.tok, self.tok_off, self.stype, self.params,
                                x, y, z, u, v, w)
