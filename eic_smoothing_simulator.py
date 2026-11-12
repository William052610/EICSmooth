#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (c) 2026 Gun Woo Lee
# EICSmooth is distributed under LGPL-3.0-only. See LICENSE and
# LICENSE-GPL-3.0.txt for the full license terms.
# This file contains C# and Java smoothing routines adapted to Python for
# EICSmooth. Adaptation notice updated 2026-11-12.
# MS-DIAL portions are derived from LGPL v3 source (Smoothing.cs).
# MZmine portions incorporate MIT-licensed code from the following files.
# SavitzkyGolayFilter.java: Copyright (c) 2004-2026 The mzmine Development Team
# SavitzkyGolaySmoothing.java: Copyright (c) 2004-2022 The MZmine Development Team
# LoessSmoothing.java: Copyright (c) 2004-2024 The mzmine Development Team
# See LICENSE-MZmine-MIT.txt for the original permission notice.
# The adapted LoessInterpolator routine includes software developed at
# The Apache Software Foundation (http://www.apache.org/).
# See LICENSE-Apache-2.0.txt, NOTICE, and THIRD_PARTY_NOTICES.md.

"""EICSmooth: LC-MS EIC smoothing comparison and morphology analysis.

This module can be launched as a desktop GUI, installed as a console command,
or imported for programmatic access to the smoothing functions.

Extract an MS1 EIC from mzML/mzXML and reproduce the smoothing implementations
used by MS-DIAL (LWMA, TBLWMA, SMA, Savitzky-Golay, binomial) and MZmine
(Savitzky-Golay, LOESS) on the same chromatogram.

Dependencies
------------
    python -m pip install -r requirements.txt

Upstream algorithm references
-----------------------------
MS-DIAL:
    https://github.com/systemsomicslab/MsdialWorkbench/blob/512f479ba2200af168d35e94bb79d1b56028b752/src/Common/CommonStandard/Algorithm/ChromSmoothing/Smoothing.cs
MZmine SG and LOESS:
    https://github.com/mzmine/mzmine/tree/51ff34834a18db599e1898b4dcae3ca5342dd51a/mzmine-community/src/main/java/io/github/mzmine/modules/dataprocessing/featdet_smoothing
Apache Commons Math 4 legacy LOESS reference for the Python port
(MZmine uses Commons Math 3.6.1):
    https://github.com/apache/commons-math/blob/f2fdc3968bdb7627029d84867398fa0af3cf0c84/commons-math-legacy/src/main/java/org/apache/commons/math4/legacy/analysis/interpolation/LoessInterpolator.java

Run
---
    python eic_smoothing_simulator.py
    python eic_smoothing_simulator.py --self-test

Peak-morphology workflow
------------------------
The raw EIC is used once to select a local maximum inside the user-defined
expected-RT search window. Candidates whose RT differs by no more than one
median scan interval are resolved by topographic prominence. The resulting
left and right peak bases are then held fixed for the raw and every smoothed
trace. This prevents sub-scan RT rounding or smoothing-specific boundary
selection from changing which chromatographic feature is quantified.

Five morphology outcomes are calculated from a trace-specific linear baseline
inside the common raw-EIC-derived boundaries: baseline-corrected peak height,
signed baseline-corrected area, FWHM, apex RT, and an explicitly defined
ideal-slope score. Raw-relative outcomes are reported as log2 height retention,
log2 area retention, log2 FWHM ratio, apex RT shift in seconds, and ideal-slope
change. No local noise region or signal-to-noise threshold is required.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import traceback
from math import comb

import numpy as np
from matplotlib.figure import Figure

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    TK_IMPORT_ERROR = None
except (ImportError, OSError) as error:
    tk = None
    filedialog = messagebox = ttk = None
    FigureCanvasTkAgg = NavigationToolbar2Tk = None
    TK_IMPORT_ERROR = error

try:
    from pyteomics import mzml, mzxml
    PYTEOMICS_IMPORT_ERROR = None
except Exception as error:
    mzml = None
    mzxml = None
    PYTEOMICS_IMPORT_ERROR = error


METHODS = ("LWMA", "TBLWMA", "SMA", "MSDIAL_SG", "BF", "MZMINE_SG", "LOESS")
METHOD_LABELS = {
    "LWMA": "MS-DIAL LWMA",
    "TBLWMA": "MS-DIAL TBLWMA",
    "SMA": "MS-DIAL SMA",
    "MSDIAL_SG": "MS-DIAL SG",
    "BF": "MS-DIAL BF",
    "MZMINE_SG": "MZmine SG",
    "LOESS": "MZmine LOESS",
}
__author__ = "Gun Woo Lee"
__version__ = "1.0"

MZMINE_SG_WIDTHS = (5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25)
CANDIDATE_RT_TIE_TOLERANCE_SCANS = 1.0

MORPHOLOGY_COLUMNS = (
    "apex_rt_min",
    "apex_rt_shift_sec",
    "peak_height_intensity",
    "log2_height_retention",
    "peak_area_intensity_sec",
    "log2_area_retention",
    "fwhm_sec",
    "log2_fwhm_ratio",
    "ideal_slope_score",
    "ideal_slope_change",
)

METRIC_CSV_FIELDS = (
    "sample_label",
    "source_file",
    "target_mz",
    "ppm",
    "trace_key",
    "trace_label",
    "software",
    "smoothing_parameter",
    "expected_rt_min",
    "reference_search_half_window_min",
    "reference_selection",
    "boundary_source",
    "search_rt_start_min",
    "search_rt_end_min",
    "analysis_rt_start_min",
    "analysis_rt_end_min",
    "reference_apex_rt_min",
    "candidate_rt_tie_tolerance_scans",
    "apex_rt_min",
    "apex_rt_shift_sec",
    "apex_intensity_raw",
    "apex_is_local_max",
    "left_base_rt_min",
    "right_base_rt_min",
    "peak_height_intensity",
    "log2_height_retention",
    "peak_area_intensity_sec",
    "log2_area_retention",
    "fwhm_sec",
    "log2_fwhm_ratio",
    "ideal_slope_score",
    "ideal_slope_change",
)

METHOD_HELP = """Software-native smoothing definitions

MS-DIAL LWMA
  Smoothing level L uses symmetric weights L-|j|+1. For L=2 the weights are (1, 2, 3, 2, 1)/9. Outside the boundary, the current point is repeated.

MS-DIAL TBLWMA
  Time-Based Linear Weighted Moving Average. It estimates the mean RT step te=(max RT-min RT)/(n-1), uses the open interval t±te(L+1), and weights an observed scan by L+1-|delta RT|/te. Only actual scans are normalized; no artificial edge padding is added.

MS-DIAL SMA
  Smoothing level L uses indices -L...+L (2L+1 points) with equal weights. Outside the data boundary, MS-DIAL repeats the current point.

MS-DIAL SG
  Smoothing level L uses 2L+1 scan positions and a cubic polynomial (four Vandermonde columns: x^0...x^3). The center row of the resulting hat matrix supplies the convolution coefficients. Outside the data boundary, MS-DIAL repeats the current point. L must be at least 2 for the cubic fit.

MS-DIAL BF
  Binomial coefficients C(2L,k), normalized across 2L+1 positions. Outside the boundary, the current point is repeated.

MZmine SG
  RT width is the full number of scans. MZmine uses fixed quadratic/cubic SG coefficients. The filter-data intersection is used at the edges without renormalization; negative results and positions originally equal to zero are set to zero. MZmine's selectable widths are 5,7,...,25 (default 5).

MZmine LOESS
  RT width is converted to bandwidth=width/n. MZmine calls Apache Commons Math LoessInterpolator with local-linear tricube weighting and 0 robustness iterations. Positions originally equal to zero remain zero (default width 5).
"""


def ppm_window(mz_value: float, ppm: float) -> tuple[float, float]:
    delta = mz_value * ppm * 1e-6
    return mz_value - delta, mz_value + delta


def parse_iso8601_duration_to_minutes(text: str) -> float | None:
    if not text or not isinstance(text, str):
        return None
    match = re.match(
        r"^P(?:T(?:(?P<h>\d+(?:\.\d+)?)H)?(?:(?P<m>\d+(?:\.\d+)?)M)?"
        r"(?:(?P<s>\d+(?:\.\d+)?)S)?)$",
        text.strip().upper(),
    )
    if not match:
        return None
    hours = float(match.group("h") or 0.0)
    minutes = float(match.group("m") or 0.0)
    seconds = float(match.group("s") or 0.0)
    return hours * 60.0 + minutes + seconds / 60.0


def extract_rt_minutes(spec: dict) -> float | None:
    try:
        scans = spec.get("scanList", {}).get("scan", [])
        if scans:
            scan = scans[0]
            value = scan.get("scan start time")
            if value is not None:
                unit = str(scan.get("unitName", "")).lower()
                return float(value) / 60.0 if "second" in unit or unit == "s" else float(value)
    except Exception:
        pass

    try:
        value = spec.get("scan start time")
        if value is not None:
            unit = str(spec.get("unitName", "")).lower()
            return float(value) / 60.0 if "second" in unit or unit == "s" else float(value)
    except Exception:
        pass

    try:
        value = spec.get("retentionTime")
        if isinstance(value, (int, float)):
            value = float(value)
            return value / 60.0 if value > 100 else value
        return parse_iso8601_duration_to_minutes(str(value))
    except Exception:
        return None


def get_ms_level(spec: dict) -> int | None:
    for key in ("ms level", "msLevel"):
        value = spec.get(key)
        if value is not None:
            try:
                return int(value)
            except (TypeError, ValueError):
                return None
    return None


def read_spectra(file_path: str):
    if mzml is None or mzxml is None:
        raise ImportError(
            "Pyteomics or one of its dependencies could not be imported.\n\n"
            f"Original error: {PYTEOMICS_IMPORT_ERROR}\n\n"
            "Install compatible dependencies with the same Python interpreter:\n"
            f'"{sys.executable}" -m pip install -r requirements.txt'
        ) from PYTEOMICS_IMPORT_ERROR
    lower = file_path.lower()
    if lower.endswith(".mzml"):
        return mzml.MzML(file_path, use_index=False)
    if lower.endswith(".mzxml"):
        return mzxml.MzXML(file_path, use_index=False)
    raise ValueError("Unsupported file format. Select an mzML or mzXML file.")


def extract_eic(
    file_path: str,
    target_mz: float,
    ppm: float = 10.0,
    rt_start: float | None = None,
    rt_end: float | None = None,
    status_callback=None,
) -> tuple[np.ndarray, np.ndarray]:
    mz_low, mz_high = ppm_window(target_mz, ppm)
    rts: list[float] = []
    intensities: list[float] = []

    with read_spectra(file_path) as reader:
        for scan_number, spec in enumerate(reader, start=1):
            ms_level = get_ms_level(spec)
            if ms_level is not None and ms_level != 1:
                continue
            rt = extract_rt_minutes(spec)
            if rt is None or (rt_start is not None and rt < rt_start):
                continue
            if rt_end is not None and rt > rt_end:
                continue

            mz_array = spec.get("m/z array")
            intensity_array = spec.get("intensity array")
            if mz_array is None or intensity_array is None:
                continue
            mz_values = np.asarray(mz_array, dtype=float)
            intensity_values = np.asarray(intensity_array, dtype=float)
            mask = (mz_values >= mz_low) & (mz_values <= mz_high)
            eic_intensity = float(np.sum(intensity_values[mask])) if np.any(mask) else 0.0
            rts.append(float(rt))
            intensities.append(eic_intensity)

            if status_callback and scan_number % 300 == 0:
                status_callback(f"Reading... {scan_number:,} scans")

    if not rts:
        return np.array([], dtype=float), np.array([], dtype=float)
    rt_array = np.asarray(rts, dtype=float)
    intensity_array = np.asarray(intensities, dtype=float)
    order = np.argsort(rt_array, kind="stable")
    return rt_array[order], intensity_array[order]


def _current_value_padded(y: np.ndarray, index: int, offset: int) -> float:
    neighbor = index + offset
    return float(y[neighbor]) if 0 <= neighbor < y.size else float(y[index])




def smooth_msdial_lwma(y: np.ndarray, smoothing_level: int) -> np.ndarray:
    """Port of MS-DIAL ``LinearWeightedMovingAverage``."""
    level = max(int(smoothing_level), 0)
    denominator = float((level + 1) ** 2)
    output = np.empty(y.size, dtype=float)
    for i in range(y.size):
        total = 0.0
        for j in range(-level, level + 1):
            weight = level - abs(j) + 1
            total += _current_value_padded(y, i, j) * weight
        output[i] = total / denominator
    return output


def smooth_msdial_tblwma(
    x: np.ndarray,
    y: np.ndarray,
    smoothing_level: int,
) -> np.ndarray:
    """Port of MS-DIAL ``TimeBasedLinearWeightedMovingAverage``."""
    if y.size <= 1:
        return y.astype(float, copy=True)
    if x.size != y.size:
        raise ValueError("TBLWMA retention-time and intensity arrays must have the same length.")
    if np.any(np.diff(x) < 0):
        raise ValueError("TBLWMA retention times must be sorted in ascending order.")

    mean_time_step = (float(np.max(x)) - float(np.min(x))) / (y.size - 1)
    if mean_time_step <= 0 or not np.isfinite(mean_time_step):
        raise ValueError("TBLWMA requires at least two distinct retention times.")

    level = max(int(smoothing_level), 0)
    distance_limit = level + 1
    output = np.empty(y.size, dtype=float)
    left_index = 0
    for i, time in enumerate(x):
        low = time - mean_time_step * distance_limit
        high = time + mean_time_step * distance_limit
        while left_index < y.size and x[left_index] <= low:
            left_index += 1

        weighted_intensity = 0.0
        total_weight = 0.0
        for k in range(left_index, y.size):
            if x[k] >= high:
                break
            weight = distance_limit - abs(time - x[k]) / mean_time_step
            weighted_intensity += float(y[k]) * weight
            total_weight += weight
        output[i] = weighted_intensity / total_weight
    return output

def smooth_msdial_sma(y: np.ndarray, smoothing_level: int) -> np.ndarray:
    """Port of MS-DIAL ``SimpleMovingAverage`` including edge behavior."""
    level = max(int(smoothing_level), 0)
    denominator = 2 * level + 1
    return np.asarray(
        [
            sum(_current_value_padded(y, i, j) for j in range(-level, level + 1))
            / denominator
            for i in range(y.size)
        ],
        dtype=float,
    )


def smooth_msdial_binomial(y: np.ndarray, smoothing_level: int) -> np.ndarray:
    """Port of MS-DIAL ``BinomialFilter`` including edge behavior."""
    level = max(int(smoothing_level), 0)
    coefficients = np.asarray([comb(2 * level, i) for i in range(2 * level + 1)], dtype=float)
    coefficients /= np.sum(coefficients)
    output = np.empty(y.size, dtype=float)
    for i in range(y.size):
        output[i] = sum(
            _current_value_padded(y, i, j) * coefficients[j + level]
            for j in range(-level, level + 1)
        )
    return output


def _msdial_sg_coefficients(smoothing_level: int) -> np.ndarray:
    """Reproduce the center row of MS-DIAL's cubic SG hat matrix."""
    level = int(smoothing_level)
    if level < 2:
        raise ValueError("The MS-DIAL SG smoothing level must be at least 2.")
    offsets = np.arange(-level, level + 1, dtype=float)
    vandermonde = np.column_stack([offsets**degree for degree in range(4)])
    gram_inverse = np.linalg.inv(vandermonde.T @ vandermonde)
    hat_matrix = vandermonde @ gram_inverse @ vandermonde.T
    return np.asarray(hat_matrix[level], dtype=float)


def smooth_msdial_sg(y: np.ndarray, smoothing_level: int) -> np.ndarray:
    """Port of MS-DIAL ``SavitxkyGolayFilter`` including edge behavior."""
    level = int(smoothing_level)
    coefficients = _msdial_sg_coefficients(level)
    output = np.empty(y.size, dtype=float)
    for i in range(y.size):
        output[i] = sum(
            _current_value_padded(y, i, j) * coefficients[j + level]
            for j in range(-level, level + 1)
        )
    return output


def _mzmine_sg_weights(width: int) -> np.ndarray:
    """Port of MZmine ``SavitzkyGolayFilter.getNormalizedWeights``."""
    width = int(width)
    if width == 0:
        return np.ones(1, dtype=float)
    if width % 2 == 0:
        width += 1
    half_width = (width - 1) // 2
    offsets = np.arange(-half_width, half_width + 1, dtype=float)
    weights = 3.0 * half_width * (half_width + 1) - 1.0 - 5.0 * offsets**2
    return weights / np.sum(weights)


def smooth_mzmine_sg(y: np.ndarray, rt_width: int) -> np.ndarray:
    """Port of MZmine ``SavitzkyGolaySmoothing.smooth``."""
    weights = _mzmine_sg_weights(rt_width)
    half_width = (weights.size - 1) // 2
    output = np.zeros(y.size, dtype=float)
    for i in range(y.size):
        k = i - half_width
        start_j = max(0, -k)
        end_j = min(weights.size, y.size - k)
        for j in range(start_j, end_j):
            output[i] += float(y[k + j]) * weights[j]
        if output[i] < 0.0:
            output[i] = 0.0
        if float(y[i]) == 0.0:
            output[i] = 0.0
    return output


def _tricube(value: float) -> float:
    absolute = abs(value)
    if absolute >= 1.0:
        return 0.0
    remainder = 1.0 - absolute**3
    return remainder**3


def _apache_loess_zero_robustness(
    x: np.ndarray,
    y: np.ndarray,
    bandwidth: float,
    accuracy: float = 1e-12,
) -> np.ndarray:
    """Port of Apache Commons Math LoessInterpolator with 0 robust iterations."""
    if x.size != y.size:
        raise ValueError("LOESS retention-time and intensity arrays must have the same length.")
    if x.size == 0:
        return np.array([], dtype=float)
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("LOESS input must not contain NaN or infinite values.")
    if np.any(np.diff(x) <= 0):
        raise ValueError(
            "As required by MZmine LOESS, retention times must be strictly "
            "increasing with no duplicates."
        )
    if x.size <= 2:
        return y.astype(float, copy=True)

    bandwidth_in_points = int(float(bandwidth) * x.size)
    if bandwidth_in_points < 2:
        raise ValueError("The MZmine LOESS RT width is too small; at least 2 scans are required.")

    output = np.empty(y.size, dtype=float)
    left = 0
    right = bandwidth_in_points - 1
    for i in range(y.size):
        if i > 0:
            next_right = right + 1
            while next_right < y.size and x[next_right] - x[i] < x[i] - x[left]:
                left += 1
                right = next_right
                next_right += 1

        edge = left if x[i] - x[left] > x[right] - x[i] else right
        denominator = abs(1.0 / (x[edge] - x[i]))
        sum_weights = 0.0
        sum_x = 0.0
        sum_x_squared = 0.0
        sum_y = 0.0
        sum_xy = 0.0
        for k in range(left, right + 1):
            distance = x[i] - x[k] if k < i else x[k] - x[i]
            weight = _tricube(distance * denominator)
            x_weighted = x[k] * weight
            sum_weights += weight
            sum_x += x_weighted
            sum_x_squared += x[k] * x_weighted
            sum_y += y[k] * weight
            sum_xy += y[k] * x_weighted

        mean_x = sum_x / sum_weights
        mean_y = sum_y / sum_weights
        mean_xy = sum_xy / sum_weights
        mean_x_squared = sum_x_squared / sum_weights
        variance_x = mean_x_squared - mean_x * mean_x
        beta = 0.0 if np.sqrt(abs(variance_x)) < accuracy else (
            (mean_xy - mean_x * mean_y) / variance_x
        )
        alpha = mean_y - beta * mean_x
        output[i] = beta * x[i] + alpha
    return output


def smooth_mzmine_loess(x: np.ndarray, y: np.ndarray, rt_width: int) -> np.ndarray:
    """Port of MZmine ``LoessSmoothing.smoothRt``."""
    if y.size <= 2:
        return y.astype(float, copy=True)
    bandwidth = min(float(rt_width) / y.size, 1.0)
    output = _apache_loess_zero_robustness(x, y, bandwidth)
    output[np.asarray(y) <= 0.0] = 0.0
    return output


def apply_smoothing(
    y: np.ndarray,
    method: str,
    msdial_level: int = 3,
    mzmine_sg_width: int = 5,
    mzmine_loess_width: int = 5,
    x: np.ndarray | None = None,
) -> tuple[np.ndarray, dict]:
    method = str(method).strip().upper()
    method = {
        "SG": "MZMINE_SG",  # backward-compatible alias used by the prior version
        "SG_MZMINE": "MZMINE_SG",
        "SG_MSDIAL": "MSDIAL_SG",
    }.get(method, method)
    y = np.asarray(y, dtype=float)
    if x is None:
        x = np.arange(y.size, dtype=float)
    else:
        x = np.asarray(x, dtype=float)

    meta = {"method": method, "label": METHOD_LABELS.get(method, method)}
    if y.size == 0:
        return y.copy(), meta
    if method in ("NONE", "RAW"):
        meta.update(source="raw", display_parameter="raw")
        return y.copy(), meta

    if method in ("LWMA", "TBLWMA", "SMA", "MSDIAL_SG", "BF"):
        level = max(int(msdial_level), 0)
        meta.update(
            source="MS-DIAL",
            smoothing_level=level,
            nominal_points=2 * level + 1,
            display_parameter=f"L={level} (2L+1={2 * level + 1})",
        )
        if method == "LWMA":
            return smooth_msdial_lwma(y, level), meta
        if method == "TBLWMA":
            return smooth_msdial_tblwma(x, y, level), meta
        if method == "SMA":
            return smooth_msdial_sma(y, level), meta
        if method == "MSDIAL_SG":
            meta["polynomial_order"] = 3
            return smooth_msdial_sg(y, level), meta
        return smooth_msdial_binomial(y, level), meta

    if method == "MZMINE_SG":
        width = int(mzmine_sg_width)
        meta.update(
            source="MZmine", rt_width_scans=width, display_parameter=f"width={width} scans"
        )
        effective_width = width + 1 if width % 2 == 0 else width
        meta["effective_width"] = effective_width
        meta["polynomial_order"] = "quadratic/cubic fixed"
        meta["display_parameter"] = f"width={effective_width} scans"
        return smooth_mzmine_sg(y, width), meta
    if method == "LOESS":
        width = int(mzmine_loess_width)
        meta.update(
            source="MZmine", rt_width_scans=width, display_parameter=f"width={width} scans"
        )
        meta["bandwidth_fraction"] = min(float(width) / y.size, 1.0)
        meta["robustness_iterations"] = 0
        return smooth_mzmine_loess(x, y, width), meta
    raise ValueError(f"Unsupported smoothing method: {method}")


def _find_local_maxima(values: np.ndarray) -> np.ndarray:
    """Return indices of non-edge local maxima, including flat-topped maxima."""
    values = np.asarray(values, dtype=float)
    if values.size < 3:
        return np.array([], dtype=int)
    left = values[1:-1] >= values[:-2]
    right = values[1:-1] >= values[2:]
    strictly_higher = (values[1:-1] > values[:-2]) | (values[1:-1] > values[2:])
    return np.flatnonzero(left & right & strictly_higher) + 1


def _interpolated_crossing(
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    target: float,
) -> float:
    if y1 == y0:
        return float((x0 + x1) / 2.0)
    fraction = (target - y0) / (y1 - y0)
    fraction = min(max(float(fraction), 0.0), 1.0)
    return float(x0 + fraction * (x1 - x0))


def _plateau_midpoint_argmax(
    values: np.ndarray,
    reference_index: int | None = None,
) -> int:
    """Return the midpoint of the most relevant exact maximum plateau.

    A single-point maximum is returned directly. If the same maximum occurs in
    several disjoint plateaus, the plateau closest to ``reference_index`` is
    used; this only resolves an otherwise exact tie and does not constrain the
    trace-specific apex search.
    """
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        raise ValueError("The peak segment used to locate the apex is empty.")
    maximum = float(np.max(values))
    maximum_indices = np.flatnonzero(values == maximum)
    runs: list[np.ndarray] = []
    run_start = 0
    for position in range(1, maximum_indices.size):
        if maximum_indices[position] != maximum_indices[position - 1] + 1:
            runs.append(maximum_indices[run_start:position])
            run_start = position
    runs.append(maximum_indices[run_start:])

    if reference_index is None:
        chosen = max(runs, key=lambda run: (run.size, -int(run[0])))
    else:
        chosen = min(
            runs,
            key=lambda run: (
                abs(float(np.mean(run)) - float(reference_index)),
                -run.size,
                int(run[0]),
            ),
        )
    return int(chosen[(chosen.size - 1) // 2])


def _log2_positive_ratio(numerator: float, denominator: float) -> float:
    if (
        not np.isfinite(numerator)
        or not np.isfinite(denominator)
        or numerator <= 0.0
        or denominator <= 0.0
    ):
        return float("nan")
    return float(np.log2(numerator / denominator))


def _validate_morphology_inputs(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size != y.size:
        raise ValueError(
            "Morphology retention-time and intensity arrays must have the same length."
        )
    if x.size < 5:
        raise ValueError("Morphology calculation requires at least 5 scans.")
    if not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("Morphology input must not contain NaN or infinite values.")
    if np.any(np.diff(x) <= 0):
        raise ValueError(
            "Morphology retention times must be strictly increasing with no duplicates."
        )
    return x, y


def _search_region_indices(
    x: np.ndarray,
    expected_rt: float | None,
    search_half_window: float | None,
) -> tuple[int, int]:
    if search_half_window is not None and search_half_window <= 0:
        raise ValueError("Reference search half-window must be greater than 0.")
    if search_half_window is not None and expected_rt is None:
        raise ValueError(
            "Expected peak RT is required when a reference search half-window is specified."
        )
    if expected_rt is not None and not (float(x[0]) <= expected_rt <= float(x[-1])):
        raise ValueError("Expected peak RT is outside the extracted EIC range.")

    if expected_rt is None or search_half_window is None:
        region_indices = np.arange(x.size, dtype=int)
    else:
        region_indices = np.flatnonzero(
            (x >= expected_rt - search_half_window)
            & (x <= expected_rt + search_half_window)
        )
    if region_indices.size < 5:
        raise ValueError("The selected reference search interval must contain at least 5 scans.")
    return int(region_indices[0]), int(region_indices[-1])


def _candidate_topographic_prominence(
    y: np.ndarray,
    candidate_index: int,
    region_start: int,
    region_end: int,
) -> tuple[float, int, int]:
    """Return raw topographic prominence and its two base indices.

    The higher of the two surrounding minima defines the contour level, as in
    the usual topographic definition of peak prominence. The minima are kept as
    fixed reference boundaries after the raw peak has been selected.
    """
    left_values = y[region_start : candidate_index + 1]
    right_values = y[candidate_index : region_end + 1]
    left_base_index = region_start + int(np.argmin(left_values))
    right_base_index = candidate_index + int(np.argmin(right_values))
    contour = max(float(y[left_base_index]), float(y[right_base_index]))
    prominence = float(y[candidate_index]) - contour
    return max(prominence, 0.0), left_base_index, right_base_index


def define_reference_peak(
    x: np.ndarray,
    raw_y: np.ndarray,
    expected_rt: float | None = None,
    search_half_window: float | None = None,
) -> dict[str, float | int | str | bool]:
    """Select one raw-EIC peak and freeze its ROI for all comparisons.

    ``expected_rt`` defines the target inside a search window. Local maxima
    separated from the closest candidate by no more than one median scan
    interval are treated as RT-equivalent; the most prominent of those
    candidates is selected. This prevents sub-scan RT rounding from selecting a
    small noise maximum while excluding more distant background peaks.
    """
    x, raw_y = _validate_morphology_inputs(x, raw_y)
    region_start, region_end = _search_region_indices(
        x, expected_rt, search_half_window
    )
    region_y = raw_y[region_start : region_end + 1]
    candidate_indices = _find_local_maxima(region_y) + region_start

    if candidate_indices.size == 0:
        apex_index = region_start + int(np.argmax(region_y))
        prominence, left_base_index, right_base_index = (
            _candidate_topographic_prominence(
                raw_y, apex_index, region_start, region_end
            )
        )
        selection = "raw_region_maximum_fallback"
        apex_is_local_max = False
    else:
        candidates: list[tuple[float, float, float, int, int, int]] = []
        for index in candidate_indices:
            prominence, left_base, right_base = _candidate_topographic_prominence(
                raw_y, int(index), region_start, region_end
            )
            distance = (
                abs(float(x[index]) - float(expected_rt))
                if expected_rt is not None
                else 0.0
            )
            candidates.append(
                (
                    prominence,
                    -distance,
                    float(raw_y[index]),
                    int(index),
                    left_base,
                    right_base,
                )
            )
        if expected_rt is None:
            selected = max(candidates)
        else:
            median_rt_step = float(np.median(np.diff(x)))
            minimum_distance = min(-candidate[1] for candidate in candidates)
            rt_equivalent = [
                candidate
                for candidate in candidates
                if -candidate[1]
                <= minimum_distance
                + CANDIDATE_RT_TIE_TOLERANCE_SCANS * median_rt_step
            ]
            selected = max(
                rt_equivalent,
                key=lambda candidate: (
                    candidate[0],
                    candidate[1],
                    candidate[2],
                ),
            )
        prominence, _, _, apex_index, left_base_index, right_base_index = selected
        selection = "expected_rt_one_scan_tie_max_raw_prominence"
        apex_is_local_max = True

    if left_base_index >= apex_index or right_base_index <= apex_index:
        raise ValueError(
            "Could not define both boundaries of the reference peak from the raw EIC. "
            "Adjust Expected peak RT or the reference search half-window."
        )

    return {
        "selection": selection,
        "search_start_index": region_start,
        "search_end_index": region_end,
        "search_rt_start_min": float(x[region_start]),
        "search_rt_end_min": float(x[region_end]),
        "expected_rt_min": (
            float(expected_rt) if expected_rt is not None else float("nan")
        ),
        "apex_index": apex_index,
        "apex_rt_min": float(x[apex_index]),
        "apex_intensity": float(raw_y[apex_index]),
        "apex_is_local_max": apex_is_local_max,
        "left_base_index": left_base_index,
        "left_base_rt_min": float(x[left_base_index]),
        "right_base_index": right_base_index,
        "right_base_rt_min": float(x[right_base_index]),
        "topographic_prominence_intensity": float(prominence),
        "candidate_rt_tie_tolerance_scans": CANDIDATE_RT_TIE_TOLERANCE_SCANS,
    }


def calculate_peak_morphology(
    x: np.ndarray,
    y: np.ndarray,
    expected_rt: float | None = None,
    half_window: float | None = None,
    reference_peak: dict[str, float | int | str | bool] | None = None,
) -> dict[str, float | int | bool | str]:
    """Calculate five absolute metrics inside raw-EIC-defined boundaries.

    If ``reference_peak`` is not supplied, it is first generated from ``y`` for
    backward-compatible single-trace use. For comparisons, callers should
    generate it once from the raw EIC with :func:`define_reference_peak` and
    pass the same object to the raw and all smoothed traces.

    A trace-specific straight baseline joins the intensities at the two fixed
    raw-derived boundary scans. The corrected trace is used for all five
    outcomes. Area is a signed trapezoidal integral in intensity-seconds.

    The ideal-slope score is an EICSmooth operational definition. Let ``I`` be
    the sum of adjacent changes in the ideal direction (up before the apex and
    down after it), and ``N`` the absolute sum in the non-ideal direction. The
    score is ``max(0, (I-N)/I)`` when ``I > 0`` and NA otherwise.

    Raw-relative values are added later by
    :func:`add_raw_relative_morphology`; they cannot be calculated from one
    trace alone.
    """
    x, y = _validate_morphology_inputs(x, y)
    if reference_peak is None:
        reference_peak = define_reference_peak(
            x,
            y,
            expected_rt=expected_rt,
            search_half_window=half_window,
        )

    required_indices = ("apex_index", "left_base_index", "right_base_index")
    if not all(key in reference_peak for key in required_indices):
        raise ValueError("The reference peak does not contain apex/base indices.")
    reference_apex_index = int(reference_peak["apex_index"])
    left_base_index = int(reference_peak["left_base_index"])
    right_base_index = int(reference_peak["right_base_index"])
    if not (
        0 <= left_base_index < reference_apex_index < right_base_index < x.size
    ):
        raise ValueError("The reference peak indices do not match the current EIC range.")

    base_x = np.asarray([x[left_base_index], x[right_base_index]], dtype=float)
    base_y = np.asarray([y[left_base_index], y[right_base_index]], dtype=float)
    segment_slice = slice(left_base_index, right_base_index + 1)
    segment_x = x[segment_slice]
    segment_y = y[segment_slice]
    baseline = np.interp(segment_x, base_x, base_y)
    corrected = segment_y - baseline

    reference_apex_segment_index = reference_apex_index - left_base_index
    apex_segment_index = _plateau_midpoint_argmax(
        corrected,
        reference_index=reference_apex_segment_index,
    )
    apex_is_local_max = bool(
        0 < apex_segment_index < corrected.size - 1
        and corrected[apex_segment_index] >= corrected[apex_segment_index - 1]
        and corrected[apex_segment_index] >= corrected[apex_segment_index + 1]
    )
    apex_index = left_base_index + apex_segment_index

    result: dict[str, float | int | bool | str] = {
        "reference_selection": str(reference_peak.get("selection", "raw_reference")),
        "boundary_source": "raw_eic_fixed",
        "search_rt_start_min": float(
            reference_peak.get("search_rt_start_min", x[left_base_index])
        ),
        "search_rt_end_min": float(
            reference_peak.get("search_rt_end_min", x[right_base_index])
        ),
        "analysis_rt_start_min": float(x[left_base_index]),
        "analysis_rt_end_min": float(x[right_base_index]),
        "expected_rt_min": float(reference_peak.get("expected_rt_min", float("nan"))),
        "reference_apex_rt_min": float(x[reference_apex_index]),
        "candidate_rt_tie_tolerance_scans": float(
            reference_peak.get(
                "candidate_rt_tie_tolerance_scans",
                CANDIDATE_RT_TIE_TOLERANCE_SCANS,
            )
        ),
        "apex_rt_min": float(x[apex_index]),
        "apex_rt_shift_sec": float("nan"),
        "apex_intensity_raw": float(y[apex_index]),
        "apex_is_local_max": apex_is_local_max,
        "left_base_rt_min": float(x[left_base_index]),
        "right_base_rt_min": float(x[right_base_index]),
        "peak_height_intensity": float("nan"),
        "log2_height_retention": float("nan"),
        "peak_area_intensity_sec": float("nan"),
        "log2_area_retention": float("nan"),
        "fwhm_sec": float("nan"),
        "log2_fwhm_ratio": float("nan"),
        "ideal_slope_score": float("nan"),
        "ideal_slope_change": float("nan"),
    }

    if left_base_index >= apex_index or right_base_index <= apex_index:
        return result

    peak_height = float(corrected[apex_segment_index])
    result["peak_height_intensity"] = peak_height
    result["peak_area_intensity_sec"] = float(
        np.sum(
            0.5
            * (corrected[:-1] + corrected[1:])
            * np.diff(segment_x)
        )
        * 60.0
    )

    epsilon = np.finfo(float).eps
    if not np.isfinite(peak_height) or peak_height <= epsilon:
        return result

    half_height = 0.5 * peak_height
    left_crossing = float("nan")
    for j in range(apex_segment_index - 1, -1, -1):
        y0 = corrected[j]
        y1 = corrected[j + 1]
        if y0 <= half_height <= y1 or y1 <= half_height <= y0:
            left_crossing = _interpolated_crossing(
                segment_x[j], segment_x[j + 1], y0, y1, half_height
            )
            break

    right_crossing = float("nan")
    for j in range(apex_segment_index, corrected.size - 1):
        y0 = corrected[j]
        y1 = corrected[j + 1]
        if y0 >= half_height >= y1 or y1 >= half_height >= y0:
            right_crossing = _interpolated_crossing(
                segment_x[j], segment_x[j + 1], y0, y1, half_height
            )
            break
    if np.isfinite(left_crossing) and np.isfinite(right_crossing) and right_crossing > left_crossing:
        result["fwhm_sec"] = float(
            (right_crossing - left_crossing) * 60.0
        )

    left_differences = np.diff(corrected[: apex_segment_index + 1])
    right_differences = -np.diff(corrected[apex_segment_index:])
    directed_differences = np.concatenate((left_differences, right_differences))
    ideal_change = float(np.sum(np.clip(directed_differences, 0.0, None)))
    nonideal_change = float(np.sum(np.clip(-directed_differences, 0.0, None)))
    if np.isfinite(ideal_change) and ideal_change > epsilon:
        result["ideal_slope_score"] = float(
            max(0.0, (ideal_change - nonideal_change) / ideal_change)
        )
    return result


def add_raw_relative_morphology(
    morphology_results: dict[str, dict[str, float | int | bool | str]],
) -> None:
    """Add the five raw-relative outcomes to morphology result dictionaries."""
    if "RAW" not in morphology_results:
        raise ValueError("RAW morphology results are required.")
    raw = morphology_results["RAW"]
    raw_height = float(raw["peak_height_intensity"])
    raw_area = float(raw["peak_area_intensity_sec"])
    raw_fwhm = float(raw["fwhm_sec"])
    raw_apex_rt = float(raw["apex_rt_min"])
    raw_slope = float(raw["ideal_slope_score"])

    for key, metrics in morphology_results.items():
        if key == "RAW":
            metrics["log2_height_retention"] = 0.0
            metrics["log2_area_retention"] = 0.0
            metrics["log2_fwhm_ratio"] = 0.0
            metrics["apex_rt_shift_sec"] = 0.0
            metrics["ideal_slope_change"] = 0.0
            continue
        metrics["log2_height_retention"] = _log2_positive_ratio(
            float(metrics["peak_height_intensity"]),
            raw_height,
        )
        metrics["log2_area_retention"] = _log2_positive_ratio(
            float(metrics["peak_area_intensity_sec"]),
            raw_area,
        )
        metrics["log2_fwhm_ratio"] = _log2_positive_ratio(
            float(metrics["fwhm_sec"]),
            raw_fwhm,
        )
        apex_rt = float(metrics["apex_rt_min"])
        metrics["apex_rt_shift_sec"] = (
            float((apex_rt - raw_apex_rt) * 60.0)
            if np.isfinite(apex_rt) and np.isfinite(raw_apex_rt)
            else float("nan")
        )
        slope = float(metrics["ideal_slope_score"])
        metrics["ideal_slope_change"] = (
            float(slope - raw_slope)
            if np.isfinite(slope) and np.isfinite(raw_slope)
            else float("nan")
        )


def analyze_morphology_traces(
    x: np.ndarray,
    raw_y: np.ndarray,
    expected_rt: float | None,
    search_half_window: float | None,
    msdial_level: int = 3,
    mzmine_sg_width: int = 5,
    mzmine_loess_width: int = 5,
    methods: tuple[str, ...] | list[str] = METHODS,
) -> tuple[
    dict[str, np.ndarray],
    dict[str, dict],
    dict[str, float | int | str | bool],
    dict[str, dict[str, float | int | bool | str]],
]:
    """Smooth one raw EIC and calculate fixed-reference morphology metrics."""
    x, raw_y = _validate_morphology_inputs(x, raw_y)
    smoothed_results: dict[str, np.ndarray] = {}
    smoothing_meta: dict[str, dict] = {}
    for method in methods:
        values, meta = apply_smoothing(
            raw_y,
            method,
            msdial_level=msdial_level,
            mzmine_sg_width=mzmine_sg_width,
            mzmine_loess_width=mzmine_loess_width,
            x=x,
        )
        smoothed_results[method] = values
        smoothing_meta[method] = meta

    reference_peak = define_reference_peak(
        x,
        raw_y,
        expected_rt=expected_rt,
        search_half_window=search_half_window,
    )
    morphology_results = {
        "RAW": calculate_peak_morphology(
            x,
            raw_y,
            reference_peak=reference_peak,
        )
    }
    for method in methods:
        morphology_results[method] = calculate_peak_morphology(
            x,
            smoothed_results[method],
            reference_peak=reference_peak,
        )
    add_raw_relative_morphology(morphology_results)
    return smoothed_results, smoothing_meta, reference_peak, morphology_results


def build_metric_rows(
    sample_label: str,
    source_file: str,
    target_mz: float,
    ppm: float,
    expected_rt: float | None,
    search_half_window: float | None,
    smoothing_meta: dict[str, dict],
    morphology_results: dict[str, dict[str, float | int | bool | str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for key, metrics in morphology_results.items():
        if key == "RAW":
            trace_label = "Raw EIC"
            software = "Raw"
            parameter = "raw"
        else:
            meta = smoothing_meta[key]
            trace_label = meta["label"]
            software = meta["source"]
            parameter = meta["display_parameter"]
        rows.append(
            {
                "sample_label": sample_label,
                "source_file": source_file,
                "target_mz": target_mz,
                "ppm": ppm,
                "trace_key": key,
                "trace_label": trace_label,
                "software": software,
                "smoothing_parameter": parameter,
                "expected_rt_min": expected_rt,
                "reference_search_half_window_min": search_half_window,
                **metrics,
            }
        )
    return rows


def write_metric_csv(path: str | os.PathLike[str], rows: list[dict[str, object]]) -> None:
    with open(path, "w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=METRIC_CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            serialized: dict[str, object] = {}
            for field in METRIC_CSV_FIELDS:
                value = row.get(field, "")
                if isinstance(value, (bool, np.bool_)):
                    serialized[field] = bool(value)
                elif isinstance(value, (float, np.floating)):
                    numeric = float(value)
                    serialized[field] = f"{numeric:.12g}" if np.isfinite(numeric) else "NA"
                else:
                    serialized[field] = value
            writer.writerow(serialized)


class EICSmoothingSimulator(tk.Tk if tk is not None else object):
    def __init__(self):
        if tk is None:
            raise ImportError(
                "The Tk GUI could not be loaded. Confirm that Python includes Tk/Tcl."
            ) from TK_IMPORT_ERROR
        super().__init__()
        self.title("EICSmooth")
        self.geometry("1480x980")
        self.minsize(1120, 820)

        self.rts: np.ndarray | None = None
        self.raw_intensities: np.ndarray | None = None
        self.smoothed_results: dict[str, np.ndarray] = {}
        self.smoothing_meta: dict[str, dict] = {}
        self.morphology_results: dict[
            str, dict[str, float | int | bool | str]
        ] = {}
        self.reference_peak: dict[str, float | int | str | bool] | None = None
        self.last_target_mz: float | None = None
        self.last_ppm: float | None = None
        self.last_label = ""

        self.file_var = tk.StringVar()
        self.mz_var = tk.StringVar()
        self.ppm_var = tk.StringVar(value="10")
        self.rt_start_var = tk.StringVar()
        self.rt_end_var = tk.StringVar()
        self.label_var = tk.StringVar()
        self.peak_rt_var = tk.StringVar()
        self.metric_half_window_var = tk.StringVar()
        self.msdial_level_var = tk.StringVar(value="3")
        self.mzmine_sg_width_var = tk.StringVar(value="5")
        self.mzmine_loess_width_var = tk.StringVar(value="5")
        self.show_raw_var = tk.BooleanVar(value=True)
        self.method_vars = {name: tk.BooleanVar(value=True) for name in METHODS}
        self.status_var = tk.StringVar(value="Ready")
        self._build_ui()

    def _build_ui(self):
        main = ttk.Frame(self, padding=10)
        main.pack(fill="both", expand=True)
        controls = ttk.LabelFrame(main, text="EIC extraction and smoothing", padding=9)
        controls.pack(fill="x")

        ttk.Label(controls, text="File").grid(row=0, column=0, sticky="w", padx=4, pady=4)
        ttk.Entry(controls, textvariable=self.file_var).grid(
            row=0, column=1, columnspan=8, sticky="ew", padx=4, pady=4
        )
        ttk.Button(controls, text="Browse", command=self.browse_file).grid(
            row=0, column=9, sticky="ew", padx=4, pady=4
        )

        field_row = ttk.Frame(controls)
        field_row.grid(row=1, column=0, columnspan=10, sticky="w", padx=4, pady=4)
        fields = (
            ("Target m/z", self.mz_var),
            ("ppm", self.ppm_var),
            ("RT start (min)", self.rt_start_var),
            ("RT end (min)", self.rt_end_var),
            ("Label", self.label_var),
        )
        for pair, (label, variable) in enumerate(fields):
            field = ttk.Frame(field_row)
            field.grid(row=0, column=pair, sticky="w", padx=(0, 12))
            ttk.Label(field, text=label).grid(row=0, column=0, sticky="w", padx=(0, 4))
            ttk.Entry(field, textvariable=variable, width=15).grid(
                row=0, column=1, sticky="w"
            )

        ttk.Label(controls, text="Expected peak RT (min)").grid(
            row=2, column=0, sticky="w", padx=4, pady=4
        )
        ttk.Entry(controls, textvariable=self.peak_rt_var, width=15).grid(
            row=2, column=1, sticky="ew", padx=4, pady=4
        )
        ttk.Label(controls, text="Reference search half-window (min)").grid(
            row=2, column=2, sticky="w", padx=4, pady=4
        )
        ttk.Entry(controls, textvariable=self.metric_half_window_var, width=15).grid(
            row=2, column=3, sticky="ew", padx=4, pady=4
        )
        ttk.Label(
            controls,
            text=(
                "Raw EIC selects a prominent peak near expected RT; "
                "its boundaries are fixed for all traces"
            ),
        ).grid(row=2, column=4, columnspan=6, sticky="w", padx=4, pady=4)

        ttk.Label(controls, text="Methods").grid(row=3, column=0, sticky="w", padx=4, pady=4)
        method_frame = ttk.Frame(controls)
        method_frame.grid(row=3, column=1, columnspan=6, sticky="w", padx=0, pady=2)
        for index, method in enumerate(METHODS):
            ttk.Checkbutton(
                method_frame, text=METHOD_LABELS[method], variable=self.method_vars[method]
            ).grid(
                row=0, column=index, padx=4
            )
        ttk.Button(method_frame, text="All", width=5, command=lambda: self._select_methods(True)).grid(
            row=0, column=len(METHODS), padx=(8, 2)
        )
        ttk.Button(method_frame, text="None", width=5, command=lambda: self._select_methods(False)).grid(
            row=0, column=len(METHODS) + 1, padx=2
        )
        ttk.Checkbutton(controls, text="Show raw", variable=self.show_raw_var).grid(
            row=3, column=7, sticky="w", padx=4
        )
        ttk.Button(controls, text="Definitions", command=lambda: messagebox.showinfo("Methods", METHOD_HELP)).grid(
            row=3, column=8, columnspan=2, sticky="ew", padx=4
        )

        ttk.Label(controls, text="MS-DIAL smoothing level").grid(
            row=4, column=0, columnspan=2, sticky="w", padx=4, pady=4
        )
        ttk.Entry(controls, textvariable=self.msdial_level_var, width=9).grid(
            row=4, column=2, sticky="ew", padx=4, pady=4
        )
        ttk.Label(controls, text="MZmine SG width").grid(
            row=4, column=3, sticky="e", padx=4
        )
        ttk.Combobox(
            controls,
            textvariable=self.mzmine_sg_width_var,
            values=MZMINE_SG_WIDTHS,
            state="readonly",
            width=7,
        ).grid(
            row=4, column=4, sticky="ew", padx=4
        )
        ttk.Label(controls, text="MZmine LOESS width").grid(
            row=4, column=5, columnspan=2, sticky="e", padx=4
        )
        ttk.Entry(controls, textvariable=self.mzmine_loess_width_var, width=9).grid(
            row=4, column=7, sticky="ew", padx=4
        )
        ttk.Button(controls, text="Extract + Compare", command=self.run_extraction).grid(
            row=5, column=5, sticky="ew", padx=4, pady=4
        )
        ttk.Button(controls, text="Reapply", command=self.reapply_smoothing).grid(
            row=5, column=6, sticky="ew", padx=4, pady=4
        )
        ttk.Button(controls, text="Save EIC CSV", command=self.save_csv).grid(
            row=5, column=7, sticky="ew", padx=4, pady=4
        )
        ttk.Button(controls, text="Save Metrics CSV", command=self.save_metrics_csv).grid(
            row=5, column=8, sticky="ew", padx=4, pady=4
        )
        ttk.Button(controls, text="Save PNG", command=self.save_plot).grid(
            row=5, column=9, sticky="ew", padx=4, pady=4
        )

        for column in range(10):
            controls.columnconfigure(column, weight=1)

        plot_frame = ttk.LabelFrame(main, text="Smoothing comparison", padding=6)
        plot_frame.pack(fill="both", expand=True, pady=(10, 0))
        self.figure = Figure(figsize=(10, 5.8), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self._reset_axes()
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(self.canvas, plot_frame, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(fill="x")

        metric_frame = ttk.LabelFrame(main, text="Peak morphology metrics", padding=6)
        metric_frame.pack(fill="x", pady=(8, 0))
        columns = (
            "method",
            "height",
            "height_log2",
            "area",
            "area_log2",
            "fwhm",
            "fwhm_log2",
            "apex_shift",
            "ideal_slope",
            "ideal_slope_change",
        )
        self.metric_table = ttk.Treeview(
            metric_frame, columns=columns, show="headings", height=8
        )
        headings = {
            "method": "Trace",
            "height": "Peak height",
            "height_log2": "log2 height retention",
            "area": "Peak area (intensity·s)",
            "area_log2": "log2 area retention",
            "fwhm": "FWHM (s)",
            "fwhm_log2": "log2 FWHM ratio",
            "apex_shift": "Apex shift (s)",
            "ideal_slope": "Ideal-slope score",
            "ideal_slope_change": "Ideal-slope change",
        }
        widths = {
            "method": 180,
            "height": 110,
            "height_log2": 135,
            "area": 155,
            "area_log2": 130,
            "fwhm": 95,
            "fwhm_log2": 125,
            "apex_shift": 105,
            "ideal_slope": 120,
            "ideal_slope_change": 135,
        }
        for column in columns:
            self.metric_table.heading(column, text=headings[column])
            self.metric_table.column(column, width=widths[column], anchor="center")
        self.metric_table.pack(fill="x", expand=True)
        ttk.Label(
            metric_frame,
            text=(
                "Peak boundaries come from the raw EIC and remain fixed. Height, signed area, "
                "FWHM, apex, and ideal slope use a trace-specific linear baseline; retention, "
                "ratio, shift, and change columns are calculated relative to RAW."
            ),
        ).pack(fill="x", pady=(4, 0))

        ttk.Label(main, textvariable=self.status_var).pack(fill="x", pady=(7, 0))

    def _reset_axes(self):
        self.ax.clear()
        self.ax.set_xlabel("Retention Time (min)")
        self.ax.set_ylabel("Intensity")
        self.ax.set_title("Extracted Ion Chromatogram")
        self.ax.grid(True, alpha=0.25)

    def _select_methods(self, selected: bool):
        for variable in self.method_vars.values():
            variable.set(selected)

    def set_status(self, text: str):
        self.status_var.set(text)
        self.update_idletasks()

    def browse_file(self):
        path = filedialog.askopenfilename(
            title="Select mzML or mzXML",
            filetypes=[
                ("Mass spectrometry files", "*.mzML *.mzXML"),
                ("mzML", "*.mzML"),
                ("mzXML", "*.mzXML"),
                ("All files", "*.*"),
            ],
        )
        if path:
            self.file_var.set(path)
            if not self.label_var.get().strip():
                self.label_var.set(os.path.splitext(os.path.basename(path))[0])

    @staticmethod
    def _float_or_none(text: str, name: str) -> float | None:
        if not text.strip():
            return None
        try:
            return float(text)
        except ValueError as error:
            raise ValueError(f"{name} must be numeric: {text}") from error

    @staticmethod
    def _integer(text: str, name: str, minimum: int = 0) -> int:
        try:
            value = int(text.strip())
        except ValueError as error:
            raise ValueError(f"{name} must be an integer: {text}") from error
        if value < minimum:
            raise ValueError(f"{name} must be at least {minimum}.")
        return value

    def _smoothing_parameters(self) -> tuple[int, int, int, list[str]]:
        msdial_level = self._integer(
            self.msdial_level_var.get(), "MS-DIAL smoothing level", 0
        )
        mzmine_sg_width = self._integer(
            self.mzmine_sg_width_var.get(), "MZmine SG RT width", 3
        )
        mzmine_loess_width = self._integer(
            self.mzmine_loess_width_var.get(), "MZmine LOESS RT width", 2
        )
        methods = [method for method in METHODS if self.method_vars[method].get()]
        if not methods and not self.show_raw_var.get():
            raise ValueError("Select at least one smoothing method or Show raw.")
        if "MSDIAL_SG" in methods and msdial_level < 2:
            raise ValueError("The MS-DIAL SG smoothing level must be at least 2.")
        if "MZMINE_SG" in methods and mzmine_sg_width not in MZMINE_SG_WIDTHS:
            raise ValueError("MZmine SG RT width must be one of 5, 7, ..., 25.")
        return msdial_level, mzmine_sg_width, mzmine_loess_width, methods

    def _metric_parameters(self) -> tuple[float | None, float | None]:
        expected_rt = self._float_or_none(self.peak_rt_var.get(), "Expected peak RT")
        half_window = self._float_or_none(
            self.metric_half_window_var.get(), "Reference search half-window"
        )
        if half_window is not None and half_window <= 0:
            raise ValueError("Reference search half-window must be greater than 0.")
        if half_window is not None and expected_rt is None:
            raise ValueError(
                "Enter Expected peak RT when using a reference search half-window."
            )
        return expected_rt, half_window

    def _extraction_parameters(self):
        path = self.file_var.get().strip()
        if not path or not os.path.isfile(path):
            raise ValueError("Select an existing mzML or mzXML file.")
        try:
            target_mz = float(self.mz_var.get().strip())
            ppm = float(self.ppm_var.get().strip())
        except ValueError as error:
            raise ValueError("Target m/z and ppm must be numeric.") from error
        if target_mz <= 0 or ppm <= 0:
            raise ValueError("Target m/z and ppm must be greater than 0.")
        rt_start = self._float_or_none(self.rt_start_var.get(), "RT start")
        rt_end = self._float_or_none(self.rt_end_var.get(), "RT end")
        if rt_start is not None and rt_end is not None and rt_start > rt_end:
            raise ValueError("RT start cannot be greater than RT end.")
        label = self.label_var.get().strip() or f"m/z {target_mz:.5f}"
        return path, target_mz, ppm, rt_start, rt_end, label

    def run_extraction(self):
        try:
            path, target_mz, ppm, rt_start, rt_end, label = self._extraction_parameters()
            self._smoothing_parameters()
            self.set_status("Extracting EIC...")
            rts, intensities = extract_eic(
                path, target_mz, ppm, rt_start, rt_end, status_callback=self.set_status
            )
            if rts.size == 0:
                raise ValueError("No matching MS1 scans were found, or EIC extraction failed.")
            self.rts = rts
            self.raw_intensities = intensities
            self.last_target_mz = target_mz
            self.last_ppm = ppm
            self.last_label = label
            self._calculate_and_plot()
        except Exception as error:
            self.set_status("Error")
            messagebox.showerror("Error", f"{error}\n\n{traceback.format_exc()}")

    def reapply_smoothing(self):
        if self.rts is None or self.raw_intensities is None:
            messagebox.showinfo("Reapply", "Run Extract + Compare first.")
            return
        try:
            self._calculate_and_plot()
        except Exception as error:
            self.set_status("Error")
            messagebox.showerror("Error", f"{error}\n\n{traceback.format_exc()}")

    def _calculate_and_plot(self):
        if self.rts is None or self.raw_intensities is None:
            return
        msdial_level, mzmine_sg_width, mzmine_loess_width, methods = (
            self._smoothing_parameters()
        )
        expected_rt, metric_half_window = self._metric_parameters()
        (
            self.smoothed_results,
            self.smoothing_meta,
            self.reference_peak,
            self.morphology_results,
        ) = analyze_morphology_traces(
            self.rts,
            self.raw_intensities,
            expected_rt=expected_rt,
            search_half_window=metric_half_window,
            msdial_level=msdial_level,
            mzmine_sg_width=mzmine_sg_width,
            mzmine_loess_width=mzmine_loess_width,
            methods=methods,
        )

        self._reset_axes()
        if self.show_raw_var.get():
            self.ax.plot(
                self.rts,
                self.raw_intensities,
                color="0.25",
                linewidth=1.0,
                alpha=0.55,
                label=f"{self.last_label} (raw)",
            )
        for method in methods:
            meta = self.smoothing_meta[method]
            self.ax.plot(
                self.rts,
                self.smoothed_results[method],
                linewidth=1.65,
                label=f"{meta['label']} ({meta['display_parameter']})",
            )
        if expected_rt is not None:
            self.ax.axvline(
                expected_rt,
                color="red",
                linewidth=1.1,
                alpha=0.8,
                label="Expected peak RT",
            )
        if self.reference_peak is not None:
            left_base_rt = float(self.reference_peak["left_base_rt_min"])
            right_base_rt = float(self.reference_peak["right_base_rt_min"])
            reference_apex_rt = float(self.reference_peak["apex_rt_min"])
            self.ax.axvspan(
                left_base_rt,
                right_base_rt,
                color="0.5",
                alpha=0.08,
                label="Fixed raw-EIC boundaries",
            )
            self.ax.axvline(
                reference_apex_rt,
                color="black",
                linestyle="--",
                linewidth=1.0,
                alpha=0.8,
                label="Raw reference apex",
            )

        title = "EIC smoothing comparison"
        if self.last_target_mz is not None and self.last_ppm is not None:
            title += f" | m/z {self.last_target_mz:.6f} ± {self.last_ppm:.2f} ppm"
        self.ax.set_title(title)
        self.ax.legend(loc="best", fontsize=9)
        self.figure.tight_layout()
        self.canvas.draw()
        self._update_metric_table(methods)

        self.set_status(
            f"Complete | points={self.rts.size:,} | methods={len(methods)} | "
            f"MS-DIAL L={msdial_level} | MZmine SG={mzmine_sg_width} | "
            f"LOESS={mzmine_loess_width} | "
            f"raw max={np.max(self.raw_intensities):.6g}"
        )

    @staticmethod
    def _format_metric(value: float | bool, digits: int = 6) -> str:
        if isinstance(value, (bool, np.bool_)):
            return "Yes" if bool(value) else "No"
        numeric = float(value)
        if not np.isfinite(numeric):
            return "NA"
        return f"{numeric:.{digits}g}"

    def _update_metric_table(self, methods: list[str]):
        for item in self.metric_table.get_children():
            self.metric_table.delete(item)
        display_order = ["RAW", *methods]
        for key in display_order:
            metrics = self.morphology_results[key]
            label = "Raw EIC" if key == "RAW" else self.smoothing_meta[key]["label"]
            self.metric_table.insert(
                "",
                "end",
                values=(
                    label,
                    self._format_metric(metrics["peak_height_intensity"]),
                    self._format_metric(metrics["log2_height_retention"]),
                    self._format_metric(metrics["peak_area_intensity_sec"]),
                    self._format_metric(metrics["log2_area_retention"]),
                    self._format_metric(metrics["fwhm_sec"]),
                    self._format_metric(metrics["log2_fwhm_ratio"]),
                    self._format_metric(metrics["apex_rt_shift_sec"]),
                    self._format_metric(metrics["ideal_slope_score"]),
                    self._format_metric(metrics["ideal_slope_change"]),
                ),
            )

    def save_csv(self):
        if self.rts is None or self.raw_intensities is None:
            messagebox.showinfo("Save CSV", "Run Extract + Compare first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save comparison CSV",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            arrays = [self.rts, self.raw_intensities]
            headers = ["rt_min", "intensity_raw"]
            for method, values in self.smoothed_results.items():
                meta = self.smoothing_meta[method]
                if meta["source"] == "MS-DIAL":
                    name = f"intensity_msdial_{method.lower()}_L{meta['smoothing_level']}"
                else:
                    name = f"intensity_mzmine_{method.lower()}_w{meta['rt_width_scans']}"
                headers.append(name)
                arrays.append(values)
            np.savetxt(
                path,
                np.column_stack(arrays),
                delimiter=",",
                header=",".join(headers),
                comments="",
                fmt="%.12g",
            )
            self.set_status(f"CSV saved: {path}")
        except Exception as error:
            messagebox.showerror("CSV save error", str(error))

    def save_metrics_csv(self):
        if self.rts is None or self.raw_intensities is None or not self.morphology_results:
            messagebox.showinfo("Save Metrics CSV", "Run Extract + Compare first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save morphology metrics CSV",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            expected_rt, metric_half_window = self._metric_parameters()
            rows = build_metric_rows(
                sample_label=self.last_label,
                source_file=os.path.basename(self.file_var.get().strip()),
                target_mz=float(self.last_target_mz),
                ppm=float(self.last_ppm),
                expected_rt=expected_rt,
                search_half_window=metric_half_window,
                smoothing_meta=self.smoothing_meta,
                morphology_results=self.morphology_results,
            )
            write_metric_csv(path, rows)
            self.set_status(f"Metrics CSV saved: {path}")
        except Exception as error:
            messagebox.showerror("Metrics CSV save error", str(error))

    def save_plot(self):
        if self.rts is None or self.raw_intensities is None:
            messagebox.showinfo("Save PNG", "Run Extract + Compare first.")
            return
        path = filedialog.asksaveasfilename(
            title="Save plot",
            defaultextension=".png",
            filetypes=[("PNG", "*.png"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            self.figure.savefig(path, dpi=600, bbox_inches="tight")
            self.set_status(f"PNG saved (600 dpi): {path}")
        except Exception as error:
            messagebox.showerror("PNG save error", str(error))


def run_self_test() -> None:
    x = np.arange(21, dtype=float)
    constant = np.full(x.size, 7.5)
    impulse = np.zeros(x.size)
    center = x.size // 2
    impulse[center] = 1.0

    for method in ("LWMA", "TBLWMA", "SMA", "MSDIAL_SG", "BF", "LOESS"):
        output, _ = apply_smoothing(
            constant,
            method,
            msdial_level=2,
            mzmine_sg_width=5,
            mzmine_loess_width=5,
            x=x,
        )
        assert output.shape == constant.shape, f"{method}: length changed"
        assert np.all(np.isfinite(output)), f"{method}: non-finite output"
        assert np.allclose(output, constant, atol=1e-9), f"{method}: constant signal changed"

    expected_lwma = np.asarray([1, 2, 3, 2, 1], dtype=float) / 9.0
    lwma, _ = apply_smoothing(impulse, "LWMA", msdial_level=2, x=x)
    tblwma, _ = apply_smoothing(impulse, "TBLWMA", msdial_level=2, x=x)
    assert np.allclose(lwma[center - 2 : center + 3], expected_lwma)
    assert np.allclose(tblwma[center - 2 : center + 3], expected_lwma)

    bf, _ = apply_smoothing(impulse, "BF", msdial_level=2, x=x)
    assert np.allclose(bf[center - 2 : center + 3], [1, 4, 6, 4, 1] / np.array(16.0))

    msdial_sg_coefficients = _msdial_sg_coefficients(2)
    expected_sg_delta = np.asarray([-3, 12, 17, 12, -3], dtype=float) / 35.0
    assert np.allclose(msdial_sg_coefficients, expected_sg_delta)

    positive_impulse = np.full(x.size, 10.0)
    positive_impulse[center] = 11.0
    msdial_sg, _ = apply_smoothing(positive_impulse, "MSDIAL_SG", msdial_level=2, x=x)
    mzmine_sg, _ = apply_smoothing(
        positive_impulse, "MZMINE_SG", mzmine_sg_width=5, x=x
    )
    assert np.allclose(msdial_sg[center - 2 : center + 3] - 10.0, expected_sg_delta)
    assert np.allclose(mzmine_sg[center - 2 : center + 3] - 10.0, expected_sg_delta)

    edge_test = np.asarray([10.0, 20.0, 30.0, 40.0, 50.0])
    msdial_edge = smooth_msdial_sg(edge_test, 2)
    expected_msdial_first = float(
        np.dot(expected_sg_delta, [edge_test[0], edge_test[0], *edge_test[:3]])
    )
    assert np.isclose(msdial_edge[0], expected_msdial_first)

    linear = 2.0 * x + 1.0
    loess, _ = apply_smoothing(linear, "LOESS", mzmine_loess_width=5, x=x)
    assert np.allclose(loess, linear, atol=1e-9)

    zero_test = np.full(x.size, 10.0)
    zero_test[center] = 0.0
    sg_zero, _ = apply_smoothing(zero_test, "MZMINE_SG", mzmine_sg_width=5, x=x)
    loess_zero, _ = apply_smoothing(
        zero_test, "LOESS", mzmine_loess_width=5, x=x
    )
    assert sg_zero[center] == 0.0 and loess_zero[center] == 0.0

    morphology_x = np.linspace(0.0, 1.0, 401)
    morphology_y = 10.0 + 100.0 * np.exp(
        -0.5 * ((morphology_x - 0.5) / 0.05) ** 2
    )
    morphology_reference = define_reference_peak(
        morphology_x,
        morphology_y,
        expected_rt=0.5,
        search_half_window=0.25,
    )
    morphology = calculate_peak_morphology(
        morphology_x,
        morphology_y,
        reference_peak=morphology_reference,
    )
    expected_fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * 0.05
    assert np.isclose(
        morphology["fwhm_sec"], expected_fwhm * 60.0, rtol=0.03
    )
    assert morphology["peak_height_intensity"] > 0.0
    assert morphology["peak_area_intensity_sec"] > 0.0
    assert np.isclose(morphology["ideal_slope_score"], 1.0)

    # Hand-calculated five-metric values for a seven-scan peak.
    metric_x = np.arange(7, dtype=float)
    metric_y = np.asarray([10, 20, 40, 80, 40, 20, 10], dtype=float)
    metric_reference = define_reference_peak(metric_x, metric_y)
    hand_metrics = calculate_peak_morphology(
        metric_x,
        metric_y,
        reference_peak=metric_reference,
    )
    assert np.isclose(hand_metrics["peak_height_intensity"], 70.0)
    assert np.isclose(hand_metrics["peak_area_intensity_sec"], 9000.0)
    assert np.isclose(hand_metrics["fwhm_sec"], 105.0)
    assert np.isclose(hand_metrics["apex_rt_min"], 3.0)
    assert np.isclose(hand_metrics["ideal_slope_score"], 1.0)

    scaled_metrics = calculate_peak_morphology(
        metric_x,
        metric_y * 0.5,
        reference_peak=metric_reference,
    )
    relative_metrics = {"RAW": hand_metrics, "SCALED": scaled_metrics}
    add_raw_relative_morphology(relative_metrics)
    assert np.isclose(scaled_metrics["log2_height_retention"], -1.0)
    assert np.isclose(scaled_metrics["log2_area_retention"], -1.0)
    assert np.isclose(scaled_metrics["log2_fwhm_ratio"], 0.0)
    assert np.isclose(scaled_metrics["apex_rt_shift_sec"], 0.0)
    assert np.isclose(scaled_metrics["ideal_slope_change"], 0.0)

    # When two candidates are RT-equivalent within one scan, a small noise
    # maximum must not displace the more prominent chromatographic peak.
    selection_x = np.linspace(0.0, 1.0, 101)
    prominent_peak = 100.0 * np.exp(
        -0.5 * ((selection_x - 0.49) / 0.010) ** 2
    )
    nearby_noise_peak = 18.0 * np.exp(
        -0.5 * ((selection_x - 0.51) / 0.004) ** 2
    )
    selection_test_y = 5.0 + prominent_peak + nearby_noise_peak
    reference_a = define_reference_peak(
        selection_x,
        selection_test_y,
        expected_rt=0.500,
        search_half_window=0.10,
    )
    reference_b = define_reference_peak(
        selection_x,
        selection_test_y,
        expected_rt=0.505,
        search_half_window=0.10,
    )
    assert reference_a["apex_index"] == reference_b["apex_index"]
    assert np.isclose(reference_a["apex_rt_min"], 0.49, atol=0.011)

    # Once the raw reference is defined, every trace must report identical
    # analysis boundaries even when smoothing changes its local apex.
    smoothed_selection_test, _ = apply_smoothing(
        selection_test_y,
        "SMA",
        msdial_level=3,
        x=selection_x,
    )
    raw_fixed = calculate_peak_morphology(
        selection_x,
        selection_test_y,
        reference_peak=reference_a,
    )
    smoothed_fixed = calculate_peak_morphology(
        selection_x,
        smoothed_selection_test,
        reference_peak=reference_a,
    )
    assert raw_fixed["left_base_rt_min"] == smoothed_fixed["left_base_rt_min"]
    assert raw_fixed["right_base_rt_min"] == smoothed_fixed["right_base_rt_min"]
    print(
        "Self-test passed: MS-DIAL LWMA/TBLWMA/SMA/SG/BF and "
        "MZmine SG/LOESS reference behaviors reproduced; prominence-based raw "
        "peak selection, fixed-boundary height/area/FWHM/apex/ideal-slope "
        "metrics, and raw-relative changes verified."
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract and compare smoothed MS1 EICs using MS-DIAL- and MZmine-derived algorithms"
    )
    parser.add_argument(
        "--self-test", action="store_true", help="run numerical checks without launching the GUI"
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args()
    if args.self_test:
        run_self_test()
        return
    app = EICSmoothingSimulator()
    app.mainloop()


if __name__ == "__main__":
    main()
