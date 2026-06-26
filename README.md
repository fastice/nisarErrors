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

## Installation

```bash
pip install git+https://github.com/fastice/nisarErrors.git@main
```

## Documentation

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
