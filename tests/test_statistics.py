import numpy as np
import pytest

from src.statistics import (
    batch_edges,
    mean_sem,
    escape_statistics,
    ratio_mean_uncertainty,
)


class TestBatchEdges:
    def test_even_split(self):
        edges = batch_edges(100, 10)
        assert list(edges) == list(range(0, 101, 10))

    def test_more_batches_than_items_is_clamped(self):
        edges = batch_edges(3, 10)
        # at most n_items batches
        assert len(edges) - 1 <= 3
        assert edges[0] == 0 and edges[-1] == 3


class TestMeanSem:
    def test_constant_values_zero_sem(self):
        mean, sem = mean_sem([5.0, 5.0, 5.0, 5.0])
        assert mean == pytest.approx(5.0)
        assert sem == pytest.approx(0.0, abs=1e-12)

    def test_known_sem(self):
        # values 1..5: mean 3, sample std sqrt(2.5), sem = sqrt(2.5)/sqrt(5)
        mean, sem = mean_sem([1, 2, 3, 4, 5])
        assert mean == pytest.approx(3.0)
        assert sem == pytest.approx(np.sqrt(2.5) / np.sqrt(5))

    def test_single_batch_sem_nan(self):
        mean, sem = mean_sem([7.0])
        assert mean == pytest.approx(7.0)
        assert np.isnan(sem)

    def test_nan_batches_dropped(self):
        mean, sem = mean_sem([2.0, np.nan, 4.0])
        assert mean == pytest.approx(3.0)
        assert not np.isnan(sem)


class TestEscapeStatistics:
    def test_uniform_leakage(self):
        # every history leaks weight 1 at 1 MeV -> leakage 1.0, energy 1e6, zero spread
        n = 100
        w = np.ones(n)
        e = np.full(n, 1e6)
        s = escape_statistics(w, e, n_source=n, n_batches=10)
        assert s['leakage'] == pytest.approx(1.0)
        assert s['leakage_sem'] == pytest.approx(0.0, abs=1e-12)
        assert s['avg_energy'] == pytest.approx(1e6)
        assert s['avg_energy_sem'] == pytest.approx(0.0, abs=1e-9)

    def test_partial_leakage_fraction(self):
        # half the histories leak -> leakage 0.5
        n = 1000
        w = np.zeros(n)
        w[::2] = 1.0
        e = np.full(n, 2e6)
        s = escape_statistics(w, e, n_source=n, n_batches=10)
        assert s['leakage'] == pytest.approx(0.5, abs=1e-9)
        assert s['avg_energy'] == pytest.approx(2e6)

    def test_weighted_average_energy(self):
        # two equally-weighted groups at 1 and 3 MeV -> mean 2 MeV
        n = 200
        w = np.ones(n)
        e = np.where(np.arange(n) < n // 2, 1e6, 3e6)
        s = escape_statistics(w, e, n_source=n, n_batches=10)
        assert s['avg_energy'] == pytest.approx(2e6, rel=1e-9)

    def test_length_mismatch_raises(self):
        with pytest.raises(ValueError):
            escape_statistics([1.0, 1.0], [1e6], n_source=2)

    def test_binned_estimator_uses_midpoints(self):
        # All escapes at 1.0 MeV; with bin edges straddling it the binned
        # estimator returns the containing bin's midpoint, not the exact value.
        n = 100
        w = np.ones(n)
        e = np.full(n, 1.0e6)
        e_bins = np.array([0.0, 0.8e6, 1.2e6, 2.0e6])  # 1 MeV -> bin [0.8,1.2], mid 1.0
        s = escape_statistics(w, e, n_source=n, n_batches=5, e_bins=e_bins)
        assert s['avg_energy'] == pytest.approx(1.0e6)

        # Shift the bin so 1 MeV lands in [0.9,1.5] (midpoint 1.2 MeV)
        e_bins2 = np.array([0.0, 0.9e6, 1.5e6, 2.0e6])
        s2 = escape_statistics(w, e, n_source=n, n_batches=5, e_bins=e_bins2)
        assert s2['avg_energy'] == pytest.approx(1.2e6)


class TestRatioMeanUncertainty:
    def test_mean_matches_weighted_average(self):
        E = [1.0, 2.0, 3.0]
        w = [1.0, 1.0, 1.0]
        s = [0.0, 0.0, 0.0]
        mean, sigma = ratio_mean_uncertainty(E, w, s)
        assert mean == pytest.approx(2.0)
        assert sigma == pytest.approx(0.0)

    def test_zero_total_weight(self):
        mean, sigma = ratio_mean_uncertainty([1, 2], [0, 0], [0, 0])
        assert mean == 0.0
        assert np.isnan(sigma)

    def test_uncertainty_propagates(self):
        # a bin far from the mean with non-zero sigma should raise sigma_E
        E = [1.0, 10.0]
        w = [100.0, 100.0]
        s = [0.0, 10.0]
        mean, sigma = ratio_mean_uncertainty(E, w, s)
        assert sigma > 0.0
