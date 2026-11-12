from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[1] / "eic_smoothing_simulator.py"
SPEC = importlib.util.spec_from_file_location("eic_smoothing_simulator", MODULE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Cannot import EICSmooth from {MODULE_PATH}")
SIMULATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SIMULATOR)


class MorphologyMetricTests(unittest.TestCase):
    def test_hand_calculated_five_metrics(self) -> None:
        x = np.arange(7, dtype=float)
        y = np.asarray([10, 20, 40, 80, 40, 20, 10], dtype=float)
        reference = SIMULATOR.define_reference_peak(x, y)
        metrics = SIMULATOR.calculate_peak_morphology(
            x,
            y,
            reference_peak=reference,
        )

        self.assertAlmostEqual(metrics["peak_height_intensity"], 70.0)
        self.assertAlmostEqual(metrics["peak_area_intensity_sec"], 9000.0)
        self.assertAlmostEqual(metrics["fwhm_sec"], 105.0)
        self.assertAlmostEqual(metrics["apex_rt_min"], 3.0)
        self.assertAlmostEqual(metrics["ideal_slope_score"], 1.0)

    def test_raw_relative_metrics(self) -> None:
        x = np.arange(7, dtype=float)
        raw_y = np.asarray([10, 20, 40, 80, 40, 20, 10], dtype=float)
        reference = SIMULATOR.define_reference_peak(x, raw_y)
        raw_metrics = SIMULATOR.calculate_peak_morphology(
            x,
            raw_y,
            reference_peak=reference,
        )
        half_metrics = SIMULATOR.calculate_peak_morphology(
            x,
            raw_y * 0.5,
            reference_peak=reference,
        )
        results = {"RAW": raw_metrics, "HALF": half_metrics}
        SIMULATOR.add_raw_relative_morphology(results)

        for field in (
            "log2_height_retention",
            "log2_area_retention",
            "log2_fwhm_ratio",
            "apex_rt_shift_sec",
            "ideal_slope_change",
        ):
            self.assertEqual(raw_metrics[field], 0.0)
        self.assertAlmostEqual(half_metrics["log2_height_retention"], -1.0)
        self.assertAlmostEqual(half_metrics["log2_area_retention"], -1.0)
        self.assertAlmostEqual(half_metrics["log2_fwhm_ratio"], 0.0)
        self.assertAlmostEqual(half_metrics["apex_rt_shift_sec"], 0.0)
        self.assertAlmostEqual(half_metrics["ideal_slope_change"], 0.0)

    def test_sloped_linear_baseline_is_removed(self) -> None:
        x = np.arange(7, dtype=float)
        shape = np.asarray([0, 10, 30, 70, 30, 10, 0], dtype=float)
        y = (5.0 + 2.0 * x) + shape
        reference = {
            "selection": "unit_test_fixed",
            "apex_index": 3,
            "left_base_index": 0,
            "right_base_index": 6,
        }
        metrics = SIMULATOR.calculate_peak_morphology(
            x,
            y,
            reference_peak=reference,
        )

        self.assertAlmostEqual(metrics["peak_height_intensity"], 70.0)
        self.assertAlmostEqual(metrics["peak_area_intensity_sec"], 9000.0)
        self.assertAlmostEqual(metrics["fwhm_sec"], 105.0)
        self.assertAlmostEqual(metrics["ideal_slope_score"], 1.0)

    def test_ideal_slope_penalizes_wrong_direction_changes(self) -> None:
        x = np.arange(7, dtype=float)
        raw_y = np.asarray([10, 20, 40, 80, 40, 20, 10], dtype=float)
        rough_y = np.asarray([10, 30, 20, 80, 40, 50, 10], dtype=float)
        reference = SIMULATOR.define_reference_peak(x, raw_y)
        raw_metrics = SIMULATOR.calculate_peak_morphology(
            x, raw_y, reference_peak=reference
        )
        rough_metrics = SIMULATOR.calculate_peak_morphology(
            x, rough_y, reference_peak=reference
        )

        self.assertAlmostEqual(raw_metrics["ideal_slope_score"], 1.0)
        self.assertLess(rough_metrics["ideal_slope_score"], 1.0)
        self.assertGreaterEqual(rough_metrics["ideal_slope_score"], 0.0)

    def test_plateau_apex_uses_middle_scan(self) -> None:
        x = np.arange(8, dtype=float)
        y = np.asarray([10, 20, 40, 80, 80, 40, 20, 10], dtype=float)
        reference = {
            "selection": "unit_test_fixed",
            "apex_index": 3,
            "left_base_index": 0,
            "right_base_index": 7,
        }
        metrics = SIMULATOR.calculate_peak_morphology(
            x,
            y,
            reference_peak=reference,
        )

        self.assertEqual(metrics["apex_rt_min"], 3.0)

    def test_fixed_boundaries_survive_smoothing(self) -> None:
        x = np.linspace(0.0, 1.0, 101)
        y = 5.0 + 100.0 * np.exp(-0.5 * ((x - 0.5) / 0.03) ** 2)
        reference = SIMULATOR.define_reference_peak(
            x,
            y,
            expected_rt=0.5,
            search_half_window=0.15,
        )
        smoothed, _ = SIMULATOR.apply_smoothing(y, "SMA", msdial_level=3, x=x)
        raw_metrics = SIMULATOR.calculate_peak_morphology(
            x, y, reference_peak=reference
        )
        smooth_metrics = SIMULATOR.calculate_peak_morphology(
            x, smoothed, reference_peak=reference
        )

        self.assertEqual(
            raw_metrics["left_base_rt_min"], smooth_metrics["left_base_rt_min"]
        )
        self.assertEqual(
            raw_metrics["right_base_rt_min"], smooth_metrics["right_base_rt_min"]
        )

    def test_nonpositive_raw_denominator_produces_na_for_ratio(self) -> None:
        self.assertTrue(np.isnan(SIMULATOR._log2_positive_ratio(1.0, 0.0)))
        self.assertTrue(np.isnan(SIMULATOR._log2_positive_ratio(-1.0, 2.0)))


if __name__ == "__main__":
    unittest.main()
