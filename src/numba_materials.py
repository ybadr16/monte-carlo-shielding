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
        u_has = []; u_emin = []; u_emax = []; u_interp = []; u_mult = []
        u_energy = []; u_e_off = [0]; u_nband = []; u_tab_off = [0]
        u_cumul = []; u_el = []; u_fis = []; u_cap = []
        in_xs = []; in_Q = []; in_mt_off = [0]; in_xs_off = [0]  # all inel MTs
        in_kind = []; in_did = []; in_nchild = []; in_pairs = []
        from .physics import watt_params_for
        nu_e = []; nu_v = []; nu_off = [0]; fissile = []; wa = []; wb = []
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
            # URR probability table (flatten cols 0/2/3/4 = cumul/el/fis/cap)
            if element not in reader._urr_cache:
                reader._load_urr(element)
            urr = reader._urr_cache.get(element)
            if urr is not None:
                eg = np.asarray(urr['energy'], float)
                tab = urr['table']
                nb = tab.shape[2]
                # multiply_smooth: bands are dimensionless factors on the smooth
                # (already 1e-24-scaled) micro XS. Otherwise they are ABSOLUTE
                # barns, so scale to cm^2 to match the 1e-24 micro convention.
                xsc = 1.0 if urr['multiply_smooth'] else 1e-24
                u_energy.extend(eg.tolist())
                for ii in range(eg.size):
                    u_cumul.extend(tab[ii, 0, :].tolist())
                    u_el.extend((tab[ii, 2, :] * xsc).tolist())
                    u_fis.extend((tab[ii, 3, :] * xsc).tolist())
                    u_cap.extend((tab[ii, 4, :] * xsc).tolist())
                u_has.append(1); u_emin.append(urr['emin']); u_emax.append(urr['emax'])
                u_interp.append(urr['interpolation']); u_mult.append(1 if urr['multiply_smooth'] else 0)
                u_nband.append(nb)
            else:
                u_has.append(0); u_emin.append(0.0); u_emax.append(0.0)
                u_interp.append(2); u_mult.append(0); u_nband.append(0)
            u_e_off.append(len(u_energy)); u_tab_off.append(len(u_cumul))
            # all neutron-emitting inelastic MTs (16/17 = n,xn; 51-90 discrete;
            # 91 continuum), micro xs on this nuclide's grid + Q. kind 0 = discrete
            # two-body, 1 = continuum secondary energy (did into the SecondaryTable);
            # nchild = extra neutrons banked ({16:1, 17:2}).
            reader._scan_inelastic_mts(element)
            gk = np.asarray(reader._fast_tables[element]["grid"], float)
            mts = [mt for mt in reader._available_inelastic_mts.get(element, [])
                   if mt in (16, 17) or 51 <= mt <= 91]
            for mt in mts:
                d = reader.get_cross_section_data(element, mt)
                if d is None:
                    continue
                xi = np.interp(gk, d["energy"], d["xs"]) * 1e-24  # micro (x N at runtime)
                in_xs.extend(xi.tolist()); in_xs_off.append(len(in_xs))
                in_Q.append(abs(d["q_value"]))
                if 51 <= mt <= 90:
                    in_kind.append(0); in_did.append(-1); in_nchild.append(0)
                else:
                    in_kind.append(1); in_did.append(len(in_pairs))
                    in_pairs.append((element, mt))
                    in_nchild.append({16: 1, 17: 2}.get(mt, 0))
            in_mt_off.append(len(in_Q))
            # nu-bar + Watt (fissile nuclides)
            reader._load_nu(element)
            nut = reader._nu_cache.get(element)
            if nut is not None:
                nu_e.extend(np.asarray(nut[0], float).tolist())
                nu_v.extend(np.asarray(nut[1], float).tolist())
                fissile.append(1)
                pa, pb = watt_params_for(element); wa.append(pa); wb.append(pb)
            else:
                fissile.append(0); wa.append(0.0); wb.append(0.0)
            nu_off.append(len(nu_e))
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
        self.urr_has = np.array(u_has, np.int64)
        self.urr_emin = np.array(u_emin); self.urr_emax = np.array(u_emax)
        self.urr_interp = np.array(u_interp, np.int64)
        self.urr_mult = np.array(u_mult, np.int64)
        self.urr_nband = np.array(u_nband, np.int64)
        self.urr_energy = np.array(u_energy) if u_energy else np.zeros(1)
        self.urr_e_off = np.array(u_e_off, np.int64)
        self.urr_tab_off = np.array(u_tab_off, np.int64)
        self.urr_cumul = np.array(u_cumul) if u_cumul else np.zeros(1)
        self.urr_el = np.array(u_el) if u_el else np.zeros(1)
        self.urr_fis = np.array(u_fis) if u_fis else np.zeros(1)
        self.urr_cap = np.array(u_cap) if u_cap else np.zeros(1)
        self.in_xs = np.array(in_xs) if in_xs else np.zeros(1)
        self.in_Q = np.array(in_Q) if in_Q else np.zeros(1)
        self.in_mt_off = np.array(in_mt_off, np.int64)   # per-nuclide MT count (cumul)
        self.in_xs_off = np.array(in_xs_off, np.int64)   # per-MT xs start in in_xs
        self.in_kind = np.array(in_kind, np.int64)       # 0 discrete two-body, 1 continuum
        self.in_nchild = np.array(in_nchild, np.int64)   # extra neutrons (n,xn)
        from .numba_secondary import SecondaryTable
        self.sec = SecondaryTable(reader.base_path, in_pairs)
        self.in_did = np.array(
            [self.sec.did_of[in_pairs[j]] if j >= 0 else -1 for j in in_did],
            np.int64)
        self.nu_e = np.array(nu_e) if nu_e else np.zeros(1)
        self.nu_v = np.array(nu_v) if nu_v else np.zeros(1)
        self.nu_off = np.array(nu_off, np.int64)
        self.fissile = np.array(fissile, np.int64)
        self.watt_a = np.array(wa); self.watt_b = np.array(wb)


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
