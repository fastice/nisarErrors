# Squint Error in Broadside-Assumed Crossing-Orbit Velocity Inversion

The two-look crossing-orbit velocity inversion derived in
[plotVerticalSensitivity.md](plotVerticalSensitivity.md) assumes each pass images at exactly
zero Doppler — i.e. that the true line-of-sight (LOS) at every pixel lies exactly broadside,
perpendicular to the platform's velocity vector. Real SAR systems deviate from this by a small
residual angle, the **squint**. This document derives how squint propagates into the velocity
solution, gives an exact closed-form result for its dominant effect, and states the correction.

It is a theory reference, not a usage guide — see the bottom of the page for the CLI invocation
of the companion plotting tool.

---

## 1. Definition and measurement of squint

Squint is the angular deviation, in the ground-projected plane, of the true LOS direction from
the broadside direction assumed by the inversion:

$$
\text{squint} = \big(\text{heading}(\hat\ell_{\text{LOS}}) - \text{heading}(\hat\ell_{\text{track}})\big)_{\pm180°} - 90°
$$

where $\hat\ell_{\text{LOS}}$ and $\hat\ell_{\text{track}}$ are the ground-projected LOS and
along-track (velocity) unit vectors, and the heading difference is wrapped to $\pm180°$ before
subtracting the nominal $90°$ broadside offset.

Squint cannot, in general, be computed from orbit state vectors alone. It depends on the
spacecraft's attitude — a data product distinct from orbit position/velocity, since nothing
about an orbit's geometry constrains which way the antenna boresight points — and, in practice,
on a residual correction empirically calibrated from the raw radar data itself (an antenna
pointed and yaw-steered exactly according to its nominal law would still leave attitude knowledge
errors that only a ground-based, data-driven calibration can resolve). Squint must therefore be
obtained either from an attitude/Doppler-centroid data product, or by direct lookup from
processor-supplied acquisition-geometry metadata, rather than re-derived from orbit mechanics.

For NISAR, this metadata is the per-pixel `losUnitVector` and `alongTrackUnitVector` fields in
each RUNW product's geolocation grid (NISAR ATBD, JPL D-95677 Rev A, §3.4–3.8) — official Level-1
processor output, already reflecting the true (squinted) acquisition geometry. Measured at scene
center for three representative ascending/descending crossing pairs (northern, central, and
southern Greenland):

| Region | squint (ascending) | squint (descending) |
|---|---|---|
| North | 1.658° | 1.494° |
| Central | 1.667° | 1.481° |
| South | 1.674° | 1.482° |

For the three Greenland tracks evaluated here, ascending squint is systematically
$\sim0.18°$ larger than descending; whether this asymmetry is a general property of the
satellite's pointing/yaw-steering behavior or specific to these tracks has not been established
from a wider sample. Within a single frame, squint varies nearly linearly with slant range (a
$\sim0.3°$–$0.4°$ spread across one swath; linear-fit residual std $\sim0.0085°$) and is nearly
constant in azimuth (std $\sim0.009°$) and height ($\sim0.013°$); across consecutive sub-frames of
the same track it drifts by only $\sim0.03°$. The behavior over a much wider heading/latitude
range than these mid-latitude examples span — in particular near a track's high-latitude turning
point, where heading changes rapidly — has not been characterized.

---

## 2. Propagation into the velocity solution

Let $\mathbf{A}_0$ be the inversion matrix (§"Geometric setup" in
[plotVerticalSensitivity.md](plotVerticalSensitivity.md)) built from the assumed broadside
headings $H_A,H_D$, and $\mathbf{A}_{\text{true}}$ the matrix built from the true,
squint-corrected headings $H_A+\text{squint}_A,\ H_D+\text{squint}_D$. For a true velocity
$\vec v_{\text{true}}$, the broadside-assumed inversion computes

$$
\vec v_{\text{computed}} = \mathbf{A}_0\,\mathbf{A}_{\text{true}}^{-1}\,\vec v_{\text{true}} = \mathbf{M}\,\vec v_{\text{true}}
$$

$\mathbf{M}$ is the identity only when squint is zero; otherwise it is the error-propagation
matrix analyzed below.

### Pure-axis error

Two reference cases avoid the direction-dependence of defining a single "% error" for an
arbitrary velocity direction: the true velocity purely along $\hat x$ or purely along $\hat y$,
giving

$$
\%\text{error}_{v_x} = (M_{00}-1)\times100, \qquad
\%\text{error}_{v_y} = (M_{11}-1)\times100
$$

### Error as a function of true flow direction

For an arbitrary true-flow orientation $\theta$, $\mathbf M\hat v(\theta)$ has both a speed error,
$(\lVert\mathbf{M}\hat v(\theta)\rVert-1)\times100$, and a direction error, the angle between
$\mathbf{M}\hat v(\theta)$ and $\hat v(\theta)$. These are generally different functions of
$\theta$ from the pure-axis cases above and from each other (§4).

---

## 3. Exact result: equal squint on both looks is a pure rotation

**Proposition.** If both looks share the same squint angle, $\text{squint}_A=\text{squint}_D=s$,
then $\mathbf{M} = R(s)$ exactly, where $R(s)$ is the $2\times2$ rotation matrix by angle $s$.

**Proof sketch.** The forward model's coefficient matrix
$N(\alpha,\beta)=\begin{pmatrix}\cos\beta&\sin\beta\\\cos(\alpha+\beta)&\sin(\alpha+\beta)\end{pmatrix}$
has rows equal to the two looks' ground-projected unit vectors, at angles $\beta$ and
$\alpha+\beta$; the inversion matrix is $\mathbf{A}=N^{-1}$. Squint rotates each look's own row by
that look's own squint angle. If both rotate by the same angle $s$, this is equivalent to right-
multiplying the whole matrix by a global rotation: $N_{\text{true}}=N_0\,R(s)$. Then
$\mathbf{M}=\mathbf{A}_0\mathbf{A}_{\text{true}}^{-1}=N_0^{-1}N_{\text{true}}=R(s)$. $\blacksquare$

This explains why, in the figures below, the dominant effect of squint on the velocity solution
is a near-constant *rotation* of the computed flow direction, close in magnitude to the mean of
the two looks' squint angles — not an arbitrary, orientation-dependent distortion. The departure
from a pure rotation is governed entirely by the *difference* $\text{squint}_A-\text{squint}_D$;
when that difference is small relative to the mean (as in the measured values above), $\mathbf M$
is close to, but not exactly, $R\big((\text{squint}_A+\text{squint}_D)/2\big)$.

---

## Figure 1 — error vs. squint magnitude, pure-axis cases

![squint error vs scale factor](images/squintError.png)

Each panel sweeps a scale factor applied to one region's measured
$(\text{squint}_A,\text{squint}_D)$ pair (1.0 reproduces the measured values exactly, preserving
the real ascending/descending asymmetry). At the measured squint:

| Region | %error $v_x$ | %error $v_y$ |
|---|---|---|
| North | +0.05% | −0.18% |
| Central | +0.08% | −0.26% |
| South | +0.05% | −0.36% |

All sub-percent; $v_y$ is consistently more sensitive than $v_x$ at these crossing geometries.

## Figure 2 — error vs. true flow direction, at the measured squint

![squint error vs flow direction](images/squintErrorByDirection.png)

Sweeping every true-flow orientation $\theta$ (0°–360°) at the measured squint: speed error
oscillates sinusoidally with $\theta$ (period $180°$), ranging from about $-0.17\%/+0.12\%$
(North) to $-0.31\%/+0.09\%$ (South). Direction error is nearly constant in $\theta$ — by the
Proposition above, this is expected, since $\mathbf M$ is close to a pure rotation — at
approximately $1.5°$–$1.7°$ for all three regions, consistent with the mean of each region's
measured ascending/descending squint.

---

## 4. The exact correction

A rotation of the *output* $(v_x,v_y)$ by $-(\text{squint}_A+\text{squint}_D)/2$ removes most of
the error when $\text{squint}_A\approx\text{squint}_D$ — for the measured values, this cuts the
maximum matrix-element error from 2.75% to 0.20%, a $\sim14\times$ reduction — but it is not the
correct general fix: it degrades as the two looks' squint values diverge, since beyond the
common-rotation part captured by the Proposition, $\mathbf M$ also has a shear component
proportional to $\text{squint}_A-\text{squint}_D$ that a single rotation cannot remove.

The exact correction instead applies each look's own measured squint to that look's own heading
**before** the inversion matrix is formed:

$$
H_A \to H_A + \text{squint}_A, \qquad H_D \to H_D + \text{squint}_D
$$

with $\alpha,\beta$ (and hence $\mathbf A$) computed from the corrected headings exactly as in the
zero-squint case. This produces $\mathbf{A}_{\text{true}}$ directly — zero residual, to the
precision of the measured squint — with no degradation as the two looks' squint values diverge,
since the divergence is exactly what gets captured rather than averaged away.

---

## 5. Application to interferometric phase vs. cross-correlation offsets

The pixel-grid position assigned to a target by zero-Doppler processing is defined purely by
geometry — the condition that the true LOS be exactly perpendicular to the true platform
velocity at the assigned time. This condition holds regardless of squint; squint is an angular
property of the LOS at that time, not a shift in which time gets assigned. Consequently, the two
sensitivity directions used to convert a measured displacement back into horizontal velocity —
the ground-projected LOS direction for range, and the ground-projected velocity direction for
azimuth — remain exactly orthogonal to each other in the data's own native geometry, independent
of squint. A system that tracks pixel displacement directly (cross-correlation/speckle tracking)
is therefore only exposed to squint through a mismatch between the *assumed* broadside heading
used in the inversion and the *true* heading the data was actually acquired with — precisely the
error quantified above, and bounded by it.

A system that forms an interferometric baseline correction has a second, independent exposure.
Converting a baseline into an equivalent range correction requires decomposing the baseline into
components parallel and perpendicular to the assumed LOS — a decomposition typically computed
from a purely two-dimensional (range, satellite-height, Earth-radius) triangle, with no
along-track or heading term at all. This implicitly assumes the LOS lies exactly in the plane
perpendicular to the ground track, an independent zero-squint assumption not shared by direct
pixel-offset tracking. Whether this second exposure is significant in a given system depends on
(a) whether the baseline itself is constructed to have a negligible along-track component
regardless of squint, and (b) the magnitude of the baseline being corrected for — a system that
has already removed the bulk of the geometric/topographic phase upstream (leaving only a small
residual baseline correction) is far less exposed than one correcting for the full physical
inter-orbit separation.

---

## Implementation notes (GrIMP `mosaic3d`)

`mosaic3d`'s zero-Doppler geometry (`common/llToImageNew.c`, `common/computeHeading.c`) assumes
exact broadside acquisition; its crossing-orbit inversion (`computeA`/`computeVxy`,
`common/initRoutines.c`) is the $\mathbf A$ of §1–§4 above. Range/azimuth offset tracking
(`make3DOffsets.c`) corresponds to §5's first case. Interferometric phase
(`computePhiZ.c`/`computePhiFlatEarth.c`) corresponds to §5's second case: the baseline-to-range
conversion uses `thetaRReZReH` (`common/initRoutines.c:918`), a pure range/Earth-radius/
satellite-height triangle with no heading term; its along-track baseline component is driven to
near zero by construction in `common/svBase.c`'s `svBaseTCN()`/`svBnBp()` (lines 271–385); and
for NISAR's ISCE flat-earth processing path the baseline being corrected is already just a small
residual orbit-error term, not the full physical separation. Combined, the current analysis
finds **no C-code change is warranted for NISAR's processing chain today** — both exposures
bound out at sub-percent speed error and the small, near-constant direction bias of §4. A
detailed, file-by-file implementation design for the exact correction (§4), including how it
would be extracted, merged across sub-frames, and threaded into the C-code chain if ever needed,
is maintained in `mosaicSource/CLAUDE.md` §"Squint (residual Doppler) sensitivity" rather than
here, since that design is specific to GrIMP's geodat/Python pipeline.

---

## Usage

```
plotSquintError [--projectDir DIR] [-o OUTPUT.png] [--outputByDirection OUTPUT2.png] [--show]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--projectDir DIR` | `/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject` | Root directory containing the example geodats and source RUNW HDF5 files. |
| `-o, --output PATH` | `squintError.png` | Output path for Figure 1 (scale sweep). |
| `--outputByDirection PATH` | `squintErrorByDirection.png` | Output path for Figure 2 (orientation sweep). |
| `--show` | off | Display interactively in addition to saving. |

## Code

[`nisarerrors/plotSquintError.py`](../nisarerrors/plotSquintError.py) — reuses `computeA`,
`headingAndIncidence`, `latLonToPS3413`, and `REGIONS` from
[`plotVerticalSensitivity.py`](../nisarerrors/plotVerticalSensitivity.py); reads squint directly
from each region's source RUNW HDF5 geolocation grid (`measureSquint()`).
