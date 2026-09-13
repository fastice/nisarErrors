# nisarErrors

Error analysis code for NISAR Ice Velocity Mapping.

**plotVerticalSensitivity** illustrates how two distinct error sources — uncompensated
vertical motion (e.g. an unmodeled tide or SMB submergence/emergence correction) and
independent range/phase measurement noise — propagate into the horizontal velocity components
$(v_x, v_y)$ produced by GrIMP's `mosaic3d` crossing-orbit InSAR-phase velocity inversion, as a
function of crossing-track geometry and pixel position. It computes real ascending/descending
track headings directly from orbital state vectors for three example crossing pairs (northern,
central, and southern Greenland) and plots both error-propagation regimes side by side for
each region.

**plotSquintError** quantifies a second-order effect: `mosaic3d` always assumes zero-Doppler
(broadside) geometry, but real acquisitions can carry a small residual squint even after
zero-Doppler processing. It computes the resulting $\%$ error in $(v_x, v_y)$ at each of the
same three real crossing-point locations, over a range of squint/PRF ratios, using a rigorous
(non-approximated) nonzero-Doppler geocoding solve.

**makeS1Reference** builds and caches a Sentinel-1 reference velocity field (the temporal mean
of the NSIDC-0725 annual Greenland mosaics, plus the inter-annual sigma and per-pixel valid-epoch
count) so that repeated calibration runs do not re-read the remote stack.

**compareToS1** compares a NISAR velocity mosaic against that reference over stable terrain and
reports whether `mosaic3d`'s formal errors describe the observed differences. Beyond the obvious
ratio it reports `kDecile` — the error ratio computed within deciles of the formal error — which
distinguishes a formal error that is merely scaled wrong (flat) from one that is missing a term
added in quadrature (declining), and a semivariogram that distinguishes short-correlation
matching noise from long-correlation ionospheric or ramp error. Those two diagnostics imply
different fixes, and a single ratio cannot tell them apart.

**velErrorStats** is the statistics library behind `compareToS1`, kept numpy/xarray-only so it
can be imported from a notebook without EarthData credentials.

## Installation

```bash
pip install git+https://github.com/fastice/nisarErrors.git@main
```

## Documentation

- [makeS1Reference](Documents/makeS1Reference.md) — what is cached and why, the grid, and the
  rioxarray/GDAL and login gotchas
- [compareToS1](Documents/compareToS1.md) — the statistics, what `kDecile` and the semivariogram
  each decide, and the point-sample vs average resolution caveat
- [plotVerticalSensitivity](Documents/plotVerticalSensitivity.md) — usage, and the
  crossing-geometry error-analysis theory (common-mode vs. independent-error propagation,
  the $\tan\gamma$ and rotation-matrix-isotropy results) behind the plots
- [plotSquintError](Documents/plotSquintError.md) — usage, and the rigorous nonzero-Doppler
  heading derivation behind the squint-error plots

## Code

- [`nisarerrors/plotVerticalSensitivity.py`](nisarerrors/plotVerticalSensitivity.py)
- [`nisarerrors/plotSquintError.py`](nisarerrors/plotSquintError.py)

## Related GrIMP documentation

This package's error analysis is grounded directly in GrIMP's `mosaic3d` C program
(`mosaicSource/Mosaic3d/make3DMosaic.c`, `mosaicSource/common/initRoutines.c`,
`mosaicSource/common/computeHeading.c`) — see `mosaicSource/Documents/mosaic3d.md` in the GrIMP
source tree for the full algorithm documentation, including the same "Error Analysis:
Crossing-Geometry Sensitivity" derivation this package's docs summarize.

## For Further Information

Please address questions to ![](https://github.com/fastice/GrIMPTools/blob/main/Email.png).
