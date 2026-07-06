"""Flat-array multi-nuclide materials for the Stage-2 Numba engine.

Compiles a list of media (each a Region with a composition of Nuclides) into
flat arrays the njit kernels traverse: per-nuclide (ragged) energy grids and
microscopic cross sections with offsets, per-nuclide A/beta and angular/
inelastic data, and a per-medium composition table (nuclide indices + number
densities). Cross sections are stored MICROSCOPIC (barns*1e-24) so a nuclide
shared across media at different number densities reuses one table.
"""
import numpy as np

from .angular_distribution import AngularDistribution

_EMAP = {"C": "C0", "Graphite": "C0", "Be": "Be9", "Al": "Al27",
         "Fe": "Fe56", "Pb": "Pb208"}


class NuclideTable:
    """Flattened per-nuclide XS / angular / inelastic data (ragged, offset-indexed)."""

    def __init__(self, reader, nuclides):
        # nuclides: list of (element, A, beta) unique across the problem
        grid, el, inl, cap, fis = [], [], [], [], []
        goff = [0]
        A = []; beta = []
        ang_eg, ang_bo, ang_mu, ang_cdf = [], [], [], []
        ang_eg_off = [0]; ang_bo_off = [0]; ang_data_off = [0]; ang_has = []
        for (element, a, b) in nuclides:
            reader._build_fast_table(element)
            tbl = reader._fast_tables[element]
            g = tbl["grid"]
            grid.extend(g.tolist())
            el.extend((tbl["el"] * 1e-24).tolist())
            inl.extend((tbl["in"] * 1e-24).tolist())
            cap.extend((tbl["cap"] * 1e-24).tolist())
            fis.extend((tbl["fis"] * 1e-24).tolist())
            goff.append(len(grid))
            A.append(a); beta.append(b)
            # angular (MT=2)
            ad = AngularDistribution(reader.base_path, _EMAP.get(element, element), mt=2)
            ad.load()
            if ad.loaded:
                eg = np.asarray(ad.energy_grid, float)
                bo = np.concatenate([np.asarray(ad.offsets, np.int64),
                                     [ad.mu_data.shape[1]]])
                ang_eg.extend(eg.tolist()); ang_eg_off.append(len(ang_eg))
                ang_bo.extend(bo.tolist()); ang_bo_off.append(len(ang_bo))
                ang_mu.extend(ad.mu_data[0].tolist())
                ang_cdf.extend(ad.mu_data[2].tolist()); ang_data_off.append(len(ang_mu))
                ang_has.append(1)
            else:
                ang_eg_off.append(len(ang_eg)); ang_bo_off.append(len(ang_bo))
                ang_data_off.append(len(ang_mu)); ang_has.append(0)
        self.grid = np.array(grid); self.goff = np.array(goff, np.int64)
        self.el = np.array(el); self.inl = np.array(inl)
        self.cap = np.array(cap); self.fis = np.array(fis)
        self.A = np.array(A); self.beta = np.array(beta)
        self.ang_eg = np.array(ang_eg) if ang_eg else np.zeros(1)
        self.ang_eg_off = np.array(ang_eg_off, np.int64)
        self.ang_bo = np.array(ang_bo, np.int64) if ang_bo else np.zeros(1, np.int64)
        self.ang_bo_off = np.array(ang_bo_off, np.int64)
        self.ang_mu = np.array(ang_mu) if ang_mu else np.zeros(1)
        self.ang_cdf = np.array(ang_cdf) if ang_cdf else np.zeros(1)
        self.ang_data_off = np.array(ang_data_off, np.int64)
        self.ang_has = np.array(ang_has, np.int64)


def compile_media(reader, mediums):
    """Return (NuclideTable, med_off, med_nuc, med_N) for a list of media whose
    Regions carry a `composition` list of Nuclide objects."""
    from .vt_calc import VelocitySampler
    uniq = {}       # element -> global nuclide index
    nuclides = []   # (element, A, beta)
    med_off = [0]; med_nuc = []; med_N = []
    for med in mediums:
        comp = med.composition
        if comp is None:
            raise ValueError("Stage-2 requires composition-based media")
        for nuc in comp:
            if nuc.element not in uniq:
                uniq[nuc.element] = len(nuclides)
                nuclides.append((nuc.element, float(nuc.atomic_weight_ratio),
                                 float(nuc.sampler.beta)))
            med_nuc.append(uniq[nuc.element])
            med_N.append(float(nuc.number_density))
        med_off.append(len(med_nuc))
    return (NuclideTable(reader, nuclides),
            np.array(med_off, np.int64), np.array(med_nuc, np.int64),
            np.array(med_N))
