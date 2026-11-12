# Third-party algorithm notices

EICSmooth includes Python adaptations of selected smoothing routines. The
links below identify the source revisions used for these adaptations. The
source file records the adaptations and their licenses.

## MS-DIAL

The LWMA, TBLWMA, SMA, Savitzky–Golay, and binomial-filter routines were
adapted from [Smoothing.cs](https://github.com/systemsomicslab/MsdialWorkbench/blob/512f479ba2200af168d35e94bb79d1b56028b752/src/Common/CommonStandard/Algorithm/ChromSmoothing/Smoothing.cs)
at commit `512f479ba2200af168d35e94bb79d1b56028b752`.
The [upstream license](https://github.com/systemsomicslab/MsdialWorkbench/blob/512f479ba2200af168d35e94bb79d1b56028b752/LGPL.txt)
is LGPL v3. The distribution includes the [LGPL v3](LICENSE) and
[GPL v3](LICENSE-GPL-3.0.txt) texts.

## MZmine

The Savitzky–Golay and LOESS routines were adapted from the MZmine
`featdet_smoothing` module at commit `51ff34834a18db599e1898b4dcae3ca5342dd51a`.
The [MIT license and copyright notices](LICENSE-MZmine-MIT.txt) are included.

| Upstream file | Copyright notice |
| --- | --- |
| [LICENSE.txt](https://github.com/mzmine/mzmine/blob/51ff34834a18db599e1898b4dcae3ca5342dd51a/LICENSE.txt) | Copyright (c) 2004-2022 The mzmine Development Team |
| [SavitzkyGolayFilter.java](https://github.com/mzmine/mzmine/blob/51ff34834a18db599e1898b4dcae3ca5342dd51a/mzmine-community/src/main/java/io/github/mzmine/modules/dataprocessing/featdet_smoothing/savitzkygolay/SavitzkyGolayFilter.java) | Copyright (c) 2004-2026 The mzmine Development Team |
| [SavitzkyGolaySmoothing.java](https://github.com/mzmine/mzmine/blob/51ff34834a18db599e1898b4dcae3ca5342dd51a/mzmine-community/src/main/java/io/github/mzmine/modules/dataprocessing/featdet_smoothing/savitzkygolay/SavitzkyGolaySmoothing.java) | Copyright (c) 2004-2022 The MZmine Development Team |
| [LoessSmoothing.java](https://github.com/mzmine/mzmine/blob/51ff34834a18db599e1898b4dcae3ca5342dd51a/mzmine-community/src/main/java/io/github/mzmine/modules/dataprocessing/featdet_smoothing/loess/LoessSmoothing.java) | Copyright (c) 2004-2024 The mzmine Development Team |

## Apache Commons Math

MZmine calls `LoessInterpolator` from [Commons Math 3.6.1](https://github.com/mzmine/mzmine/blob/51ff34834a18db599e1898b4dcae3ca5342dd51a/gradle/libs.versions.toml).
EICSmooth's Python LOESS routine follows the zero-robustness path in
[LoessInterpolator.java](https://github.com/apache/commons-math/blob/f2fdc3968bdb7627029d84867398fa0af3cf0c84/commons-math-legacy/src/main/java/org/apache/commons/math4/legacy/analysis/interpolation/LoessInterpolator.java)
from the Commons Math 4 legacy module at commit
`f2fdc3968bdb7627029d84867398fa0af3cf0c84`. The [Apache License 2.0](LICENSE-Apache-2.0.txt)
and [upstream notices](NOTICE) are included.

The five morphology metrics and their raw-relative transformations are
EICSmooth definitions documented in `docs/METRICS.md`. They do not copy or
translate a third-party morphology implementation. The term "ideal slope" is
used descriptively for monotonic rise-to-apex and fall-from-apex behavior; the
normalization used here is stated explicitly and should not be interpreted as
an exact reproduction of an MS-DIAL peak-quality score.
