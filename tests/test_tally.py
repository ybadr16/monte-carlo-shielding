import numpy as np
import pytest

from src.tally import Tally


def _result(final_weight=0.0, final_energy=0.0, absorbed_weight=0.0,
            absorbed_coords=None, region_detected=False, result="simulated"):
    return {
        "result": result,
        "absorbed_weight": absorbed_weight,
        "absorbed_coords": absorbed_coords or [],
        "final_energy": final_energy,
        "final_weight": final_weight,
        "region_detected": region_detected,
    }


class TestTally:
    def test_leakage_and_energy_statistics(self):
        tally = Tally()
        n = 100
        for _ in range(n):
            tally.merge_partial_results(_result(final_weight=1.0, final_energy=2e6))
        stats = tally.statistics(num_particles=n, n_batches=10)
        assert stats['leakage'] == pytest.approx(1.0)
        assert stats['leakage_sem'] == pytest.approx(0.0, abs=1e-12)
        assert stats['avg_energy'] == pytest.approx(2e6)

    def test_partial_leakage(self):
        tally = Tally()
        n = 100
        for i in range(n):
            leaked = 1.0 if i % 2 == 0 else 0.0
            tally.merge_partial_results(_result(final_weight=leaked, final_energy=1e6))
        stats = tally.statistics(num_particles=n, n_batches=10)
        assert stats['leakage'] == pytest.approx(0.5, abs=1e-9)

    def test_absorbed_weight_accumulates(self):
        tally = Tally()
        tally.merge_partial_results(_result(absorbed_weight=0.3))
        tally.merge_partial_results(_result(absorbed_weight=0.7))
        assert tally.results["absorbed"] == pytest.approx(1.0)

    def test_get_results_no_crash(self):
        # Regression: get_results used to reference a non-existent attribute.
        tally = Tally()
        tally.merge_partial_results(_result(final_weight=1.0, final_energy=1e6,
                                            absorbed_coords=[(0, 0, 0, 0.1)]))
        out = tally.get_results()
        assert set(out) >= {"results", "absorbed_coordinates",
                            "energy_spectrum", "leaked_weights", "region_count"}

    def test_print_summary_runs(self, capsys):
        tally = Tally()
        for _ in range(20):
            tally.merge_partial_results(_result(final_weight=1.0, final_energy=1e6))
        tally.print_summary(20, n_batches=5)
        out = capsys.readouterr().out
        assert "Leakage Fraction" in out
        assert "±" in out
