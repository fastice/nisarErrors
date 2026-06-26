# plotSquintError

Quantifies the error in `mosaic3d`'s computed $(v_x, v_y)$ caused by ignoring residual squint
(nonzero Doppler) at each of `plotVerticalSensitivity`'s three real crossing-point locations
(northern/central/southern Greenland), evaluated only at each region's real pixel position —
not swept over pixel position the way the main sensitivity plot is.

---

## Usage

```
plotSquintError [--projectDir DIR] [-o OUTPUT.png] [--show]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--projectDir DIR` | `/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject` | Root directory containing the example geodats (same `REGIONS` as `plotVerticalSensitivity`). |
| `-o, --output PATH` | `squintError.png` | Output plot path. |
| `--show` | off | Display interactively in addition to saving. |

---

## Theory

`mosaic3d` always assumes zero-Doppler (broadside) geometry when computing heading
(`common/computeHeading.c`), so it builds its solving matrix from the *assumed* zero-squint
headings: $\mathbf{A}_0 = \text{computeA}(H_A(0), H_D(0))$. The *true* measurement reflects
whatever the actual acquisition geometry was — possibly squinted:
$\mathbf{A}_{\text{true}} = \text{computeA}(H_A(f_{dop,A}), H_D(f_{dop,D}))$. For a true
velocity $\vec{v}_{\text{true}}$, the true measurements are
$(p_A,p_D) = \mathbf{A}_{\text{true}}^{-1}\vec{v}_{\text{true}}$ (flat terrain, $\mathbf{B}=0$),
and mosaic3d computes:

$$
\vec{v}_{\text{computed}} = \mathbf{A}_0\,(p_A,p_D) = \underbrace{\mathbf{A}_0\,
\mathbf{A}_{\text{true}}^{-1}}_{\mathbf{M}}\,\vec{v}_{\text{true}}
$$

Since $\mathbf{M}$ is linear, the resulting fractional error doesn't depend on
$|\vec{v}_{\text{true}}|$, only its direction. Each component's error is reported against its
own natural reference case (true velocity purely along that axis) — well-defined, no
direction assumption beyond that, no division-by-zero:

$$
\%\text{error}_{v_x} = (M_{00}-1)\times100 \quad(\vec{v}_{\text{true}}=\hat{x}), \qquad
\%\text{error}_{v_y} = (M_{11}-1)\times100 \quad(\vec{v}_{\text{true}}=\hat{y})
$$

### Computing heading at nonzero Doppler — rigorously, not by approximation

A naive small-angle estimate (antenna squint angle
$\theta=\arcsin(f_{dop}\lambda/2V_{sat})$ added directly to heading) was checked numerically
against a rigorous calculation and found wrong by a factor of ~4 — the actual heading shift is
real and roughly linear in squint, but not equal to the raw antenna squint angle. This script
therefore computes heading at a target Doppler **rigorously**: `geodatrxa.llzPtToRA()`'s
zero-Doppler Newton solve, $(\vec{T}-\vec{S})\cdot\vec{V} = 0$, is generalized to a nonzero
target,

$$
(\vec{T}-\vec{S})\cdot\vec{V} = \frac{f_{dop}\,\lambda\,R}{2}
$$

(derived from $f_{dop} = -\frac{2}{\lambda}\frac{dR}{dt} = \frac{2}{\lambda}
\frac{(\vec{T}-\vec{S})\cdot\vec{V}}{R}$), reusing only `geodatrxa`'s Doppler-agnostic
primitives (`lltoecef`, `interpPos`, `interpVel`, `rNearSLP`, `slpRg`, `prf`) — `geodatrxa.py`
itself is not modified.

The squint/PRF ratio (x-axis, $\pm2$) is applied identically to both the ascending and
descending image, each using its own real PRF — i.e. a shared fractional residual squint, not
independently-varying squints per image. The secondary top x-axis shows the **true antenna
squint angle** ($\theta=\arcsin(f_{dop}\lambda/2V_{sat})$, region-specific real PRF/$\lambda$/
$V_{sat}$) — a well-defined physical quantity, distinct from (and not simply proportional in a
known closed form to) the resulting heading shift.

### What's *not* included

Only the heading ($\mathbf{A}$-matrix) error is captured. `computeA()` does not depend on
incidence angle $\psi$ at all (confirmed against `common/initRoutines.c`), so a separate,
smaller effect (~5–6× less, confirmed numerically) on the $\psi$-dependent phase/offset-to-
velocity-units scaling is excluded — whether mosaic3d's geodat-metadata $\psi$ already
reflects the true acquisition geometry is ambiguous without knowing the NISAR processor's
internals.

## Code

[`nisarerrors/plotSquintError.py`](../nisarerrors/plotSquintError.py) — reuses `computeA`,
`headingAndIncidence`, `latLonToPS3413`, and `REGIONS` from
[`plotVerticalSensitivity.py`](../nisarerrors/plotVerticalSensitivity.py).
