# EICSmooth morphology definitions

## Common peak region and baseline

For each replicate, the unsmoothed EIC is used to select a reference peak. Its
left and right boundary scan indices, `L` and `R`, are then reused for RAW and
all smoothed traces. Let trace `m` have intensity `y(m,i)` at retention time
`t(i)`. The trace-specific straight baseline is

```text
b(m,i) = y(m,L) + [y(m,R)-y(m,L)] [t(i)-t(L)] / [t(R)-t(L)]
```

and the corrected intensity is

```text
z(m,i) = y(m,i) - b(m,i).
```

The apex `k(m)` is the maximum `z(m,i)` inside `L...R`. For a flat maximum,
EICSmooth uses the middle scan of the plateau (the earlier middle scan when the
plateau has an even number of points). Exact tied plateaus are resolved by
proximity to the raw reference apex.

## 1. Peak-height retention

```text
H(m) = z(m,k(m))
log2 height retention(m) = log2[H(m) / H(RAW)]
```

The absolute unit is intensity. The ratio is `NA` if either height is
nonpositive or nonfinite.

## 2. Peak-area retention

```text
A(m) = 60 * sum(i=L...R-1) {
  [z(m,i) + z(m,i+1)] / 2 * [t(i+1)-t(i)]
}
log2 area retention(m) = log2[A(m) / A(RAW)]
```

Retention time is supplied in minutes and the factor 60 gives an output unit of
intensity-seconds. The area is signed: negative corrected values are not
clipped. The ratio is `NA` if either area is nonpositive or nonfinite.

## 3. FWHM ratio

The half-height is `H(m)/2`. EICSmooth searches outward from the apex for the
first crossing on each side and linearly interpolates between adjacent scans.

```text
FWHM(m) = 60 * [t(right crossing) - t(left crossing)]
log2 FWHM ratio(m) = log2[FWHM(m) / FWHM(RAW)]
```

The absolute unit is seconds. FWHM is `NA` if either crossing is absent, and the
ratio is `NA` if either FWHM is nonpositive or nonfinite.

## 4. Apex RT shift

```text
apex RT shift(m) = 60 * [t(k(m)) - t(k(RAW))]
```

The unit is seconds. Positive values indicate a later apex and negative values
an earlier apex.

## 5. Ideal-slope change

For the left side, an ideal adjacent change is an increase toward the apex. For
the right side, an ideal adjacent change is a decrease away from the apex. The
directed changes are therefore

```text
d(i) = z(m,i+1) - z(m,i)        for L <= i < k(m)
d(i) = z(m,i) - z(m,i+1)        for k(m) <= i < R.
```

Let

```text
I(m) = sum max[d(i), 0]
N(m) = sum max[-d(i), 0].
```

The EICSmooth ideal-slope score is

```text
S(m) = max{0, [I(m)-N(m)] / I(m)}.
```

It ranges from 0 to 1 and is `NA` when `I(m)` is zero or nonfinite. The
raw-relative outcome is

```text
ideal-slope change(m) = S(m) - S(RAW).
```

This is an operational EICSmooth definition. It describes the same ideal
directionality concept used for chromatographic model peaks, but its stated
normalization should not be labeled as an exact MS-DIAL native score.

## RAW row and replicate handling

All five raw-relative fields are set to 0 in the RAW row. Metrics are calculated
for each replicate EIC separately. EICSmooth does not average replicate traces
before measurement; downstream summaries can report the median and range or
interquartile range for each compound-by-method combination.

## Methods-ready text

Peak morphology was quantified separately for each replicate EIC. A reference
peak and its left and right boundaries were selected from the unsmoothed EIC,
and the same boundary scan indices were applied to the unsmoothed and all
smoothed traces. For each trace, a linear local baseline connecting the two
boundary intensities was subtracted. Baseline-corrected peak height, signed
trapezoidal area, full width at half maximum (FWHM), apex retention time, and an
ideal-slope score were calculated within the common boundary. Height, area, and
FWHM were expressed as log2 ratios to the unsmoothed EIC; apex displacement was
expressed in seconds; and ideal-slope change was calculated by subtraction of
the unsmoothed score. The ideal-slope score was defined as max[0, (I-N)/I],
where I and N denote the cumulative absolute scan-to-scan intensity changes in
the ideal and nonideal directions, respectively. Values requiring nonpositive
ratio denominators or missing half-height crossings were reported as not
available.
