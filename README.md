# EICSmooth

EICSmooth is a Python GUI tool for comparing smoothing-induced changes in
LC-MS extracted-ion chromatograms (EICs). It applies seven smoothing methods
to the same MS1 EIC and quantifies peak morphology using one set of peak
boundaries defined from the unsmoothed trace.

EICSmooth was developed for the manuscript **Smoothing Alters Feature
Composition and Low-Intensity Signal Recovery in Liquid Chromatography–Mass
Spectrometry (LC–MS) Untargeted Metabolomics**.

## Included smoothing methods

- MS-DIAL: linear weighted moving average (LWMA), time-based LWMA (TBLWMA),
  simple moving average (SMA), Savitzky-Golay (SG), and binomial filter (BF)
- MZmine: Savitzky-Golay (SG) and LOESS smoothing (LS)

The source snapshots used for the Python ports are pinned in the module and
documented in `THIRD_PARTY_NOTICES.md`.

## Five morphology outcomes

All traces use the same peak boundaries selected from the raw EIC. Within
those boundaries, each trace receives its own straight baseline joining the
intensities at the two boundaries.

| Absolute value | RAW-relative value | Interpretation |
| --- | --- | --- |
| `peak_height_intensity` | `log2_height_retention` | Apex height above the local baseline |
| `peak_area_intensity_sec` | `log2_area_retention` | Signed baseline-corrected trapezoidal area |
| `fwhm_sec` | `log2_fwhm_ratio` | Width at half of the baseline-corrected height |
| `apex_rt_min` | `apex_rt_shift_sec` | Apex position and its shift from RAW |
| `ideal_slope_score` | `ideal_slope_change` | Monotonic rise-to-apex and fall-from-apex regularity |

The three ratios are reported on a log2 scale. A value of 0 means no change,
-1 means one-half of RAW, and +1 means twice RAW. The ideal-slope score ranges
from 0 to 1; 1 indicates no scan-to-scan movement in the wrong direction.

S/N and threshold-derived effective scan counts are intentionally excluded,
so the workflow does not require a user-defined noise region. Exact equations,
units, invalid-value rules, and a Methods-ready description are provided in
`docs/METRICS.md`.

## Installation on Windows

EICSmooth supports Python 3.10-3.12. A project-local virtual
environment is recommended because it keeps EICSmooth separate from packages
already installed in Anaconda `base` or another Python environment.

Open PowerShell or the PyCharm Terminal and move to the extracted EICSmooth
folder:

```powershell
cd "C:\path\to\EICSmooth"
```

Confirm that the selected Python version is supported:

```bash
python --version
```

Create the environment and install EICSmooth:

```bash
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

The last command installs EICSmooth itself and the tested runtime versions:
NumPy 1.26.4, Matplotlib 3.9.3, and Pyteomics 5.0.1 with XML support.

## GUI

Start EICSmooth with the Python executable inside the project environment:

```bash
.\.venv\Scripts\python.exe eic_smoothing_simulator.py
```

Activation is not required. In PyCharm, the same environment can be selected
as the project interpreter at `.venv\Scripts\python.exe`; after that, the Run
button can launch `eic_smoothing_simulator.py` directly.

1. Select an mzML or mzXML file.
2. Enter the target m/z, ppm tolerance, and EIC RT range.
3. Enter the expected peak RT and reference search half-window.
4. Select the smoothing methods and parameters.
5. Click **Extract + Compare**.
6. Save the EIC traces, morphology metrics, or 600-dpi PNG.

Entering the expected peak RT and search half-window is strongly recommended.
If omitted, EICSmooth searches the entire extracted EIC and may select a
different chromatographic peak.

## Synthetic example

The synthetic example does not require an mzML or mzXML file:

```bash
.\.venv\Scripts\python.exe examples\run_synthetic_example.py
```

It generates a deterministic EIC, an eight-trace morphology table, and a
600-dpi comparison figure in `examples/output/`.



## Scope and limitations

- EICSmooth compares EIC transformations; it is not a feature detector or a
  metabolite identification program.
- Peak boundaries are an operational raw-EIC reference and should be inspected
  for the compounds reported in a manuscript.
- Signed area retains values below the straight local baseline; methods that
  produce a nonpositive area are reported as `NA` for log2 area retention.
- FWHM is `NA` when either half-height crossing is absent inside the fixed
  boundaries.
- The ideal-slope score is the explicit EICSmooth definition in
  `docs/METRICS.md`, not an exact reproduction of a software-native peak-quality
  score.

## Citation

Author: Gun Woo Lee, Konkuk University  
Contact: [dlrjsdn0526@konkuk.ac.kr](mailto:dlrjsdn0526@konkuk.ac.kr)  
Repository: <https://github.com/William052610/EICSmooth>

## License

EICSmooth is distributed under the GNU Lesser General Public License v3.0.
See the [LGPL v3](LICENSE) and [GPL v3](LICENSE-GPL-3.0.txt) license texts and
the [third-party notices](THIRD_PARTY_NOTICES.md) for the adapted routines.

Copyright (c) 2026 Gun Woo Lee.
