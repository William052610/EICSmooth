#!/usr/bin/env python3
"""Generate and analyze a deterministic synthetic EIC with EICSmooth."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import eic_smoothing_simulator as eicsmooth  # noqa: E402


def synthetic_eic() -> tuple[np.ndarray, np.ndarray]:
    rt = np.linspace(4.70, 5.30, 121)
    baseline = 120.0 + 35.0 * (rt - 5.00)
    peak = 2200.0 * np.exp(-0.5 * ((rt - 5.00) / 0.035) ** 2)
    structured_noise = (
        32.0 * np.sin(2.0 * np.pi * np.arange(rt.size) / 3.0)
        + 15.0 * np.cos(2.0 * np.pi * np.arange(rt.size) / 7.0)
    )
    intensity = np.maximum(baseline + peak + structured_noise, 0.0)
    return rt, intensity


def write_eic_csv(path: Path, rt: np.ndarray, intensity: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rt_min", "intensity_raw"])
        for rt_value, intensity_value in zip(rt, intensity):
            writer.writerow([f"{rt_value:.6f}", f"{intensity_value:.12g}"])


def main() -> None:
    example_dir = Path(__file__).resolve().parent
    output_dir = example_dir / "output"
    output_dir.mkdir(exist_ok=True)

    rt, raw = synthetic_eic()
    write_eic_csv(example_dir / "synthetic_eic.csv", rt, raw)
    smoothed, meta, reference, metrics = eicsmooth.analyze_morphology_traces(
        rt,
        raw,
        expected_rt=5.00,
        search_half_window=0.20,
        msdial_level=3,
        mzmine_sg_width=5,
        mzmine_loess_width=5,
    )
    rows = eicsmooth.build_metric_rows(
        sample_label="synthetic_peak",
        source_file="synthetic_eic.csv",
        target_mz=200.0000,
        ppm=10.0,
        expected_rt=5.00,
        search_half_window=0.20,
        smoothing_meta=meta,
        morphology_results=metrics,
    )
    metric_path = output_dir / "synthetic_metrics.csv"
    eicsmooth.write_metric_csv(metric_path, rows)

    figure = Figure(figsize=(9.0, 5.5), dpi=120)
    FigureCanvasAgg(figure)
    axis = figure.add_subplot(111)
    axis.plot(rt, raw, color="0.2", linewidth=1.0, alpha=0.7, label="RAW")
    for method in eicsmooth.METHODS:
        axis.plot(rt, smoothed[method], linewidth=1.2, label=meta[method]["label"])
    axis.axvspan(
        float(reference["left_base_rt_min"]),
        float(reference["right_base_rt_min"]),
        color="0.5",
        alpha=0.08,
        label="Fixed RAW boundary",
    )
    axis.set_xlabel("Retention time (min)")
    axis.set_ylabel("Intensity")
    axis.set_title("EICSmooth synthetic example")
    axis.grid(alpha=0.2)
    axis.legend(fontsize=7, ncol=2)
    figure.tight_layout()
    figure.savefig(output_dir / "synthetic_comparison.png", dpi=600)

    print(f"Wrote {example_dir / 'synthetic_eic.csv'}")
    print(f"Wrote {metric_path}")
    print(f"Wrote {output_dir / 'synthetic_comparison.png'}")


if __name__ == "__main__":
    main()
