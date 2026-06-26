# plotVerticalSensitivity

Illustrates how two distinct error sources propagate into the horizontal velocity components
$(v_x, v_y)$ produced by `mosaic3d`'s crossing-orbit InSAR-phase velocity inversion
(`computeA`/`computeVxy` in GrIMP's `mosaicSource/common/initRoutines.c`), as a function of
crossing-track geometry and pixel position. Computes real ascending/descending track headings
from orbital state vectors for three example crossing pairs (northern, central, and southern
Greenland), and plots both error-propagation regimes side by side for each.

---

## Usage

```
plotVerticalSensitivity [--projectDir DIR] [-o OUTPUT.png] [--show]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--projectDir DIR` | `/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject` | Root directory containing the `track-N/orbit_NNNN/geodat26x16.geojson` example files. |
| `-o, --output PATH` | `verticalSensitivity.png` | Output plot path. |
| `--show` | off | Display the figure interactively in addition to saving it. |

The three example crossing pairs (track/frame numbers, relative to `--projectDir`) are fixed
in `REGIONS` at the top of `plotVerticalSensitivity.py` — edit that list directly to point at
different track/frame examples for a different project; only the root directory is a CLI
option.

| Region | Ascending geodat | Descending geodat |
|--------|-------------------|---------------------|
| Northern Greenland | `track-58/3271_37` | `track-97/3656_53` |
| Central Greenland | `track-131/3517_0001` | `track-11/3397_0000` |
| Southern Greenland | `track-45/3431_0000` | `track-126/3512_0002` |

Each region's ascending/descending heading ($H_A$, $H_D$), incidence angle, and real pixel
position ($\gamma$, marked on the plot) are computed at runtime directly from orbital state
vectors in the geodat files and printed to stdout.

---

## Theory: crossing-geometry error sensitivity

Both panels solve the same underlying 2×2 inversion `mosaic3d` uses to combine an ascending and
a descending line-of-sight (LOS) measurement into horizontal velocity:

$$
\begin{pmatrix} v_x \\ v_y \end{pmatrix} = \mathbf{A}
\begin{pmatrix} p_A \\ p_D \end{pmatrix}, \qquad
\mathbf{A} = \frac{1}{\sin^2\!\alpha}
\begin{pmatrix}
\cos\beta - \cos\alpha\cos(\alpha+\beta) & \cos(\alpha+\beta) - \cos\alpha\cos\beta \\
\sin\beta - \cos\alpha\sin(\alpha+\beta) & \sin(\alpha+\beta) - \cos\alpha\sin\beta
\end{pmatrix}
$$

(flat-terrain approximation, $\mathbf{B}=0$ — see `mosaic3d.md`'s "Matrix B" for the neglected
DEM-slope term), where $\alpha = H_A - H_D$ is the heading difference between the two tracks,
$\beta = \phi - H_A$ is the pixel's azimuth angle (in map coordinates) relative to the
ascending heading, and $H_A, H_D$ are each track's **cross-track** heading (not its
flight-direction heading — these differ by a curved-Earth, incidence-angle-dependent amount,
confirmed ~6–16° for real NISAR geometry; see `headingAndIncidence()`'s docstring in the code).

Define $H_{\text{mean}} = (H_A+H_D)/2$ and $\gamma = \phi - H_{\text{mean}}$ (pixel azimuth
relative to the *mean* heading, plotted on the x-axis) — the two error sources depend on
$\alpha$ and $\gamma$ in qualitatively different ways:

### 1. Common-mode errors (right panel: uncompensated $v_z^{SMB}$)

The same physical vertical rate $v_z^{SMB}$ contributes to both LOS measurements, scaled by
each image's own incidence angle: $\Delta p_A = v_z^{SMB}\cot\psi_A$,
$\Delta p_D = v_z^{SMB}\cot\psi_D$ (not simply $v_z^{SMB}\times1$, which only holds at
$\psi=45°$). For equal contributions $\delta$ to both measurements, $\mathbf{A}\cdot(\delta,
\delta)^T$ reduces via sum-to-product identities to a clean closed form:

$$
\Delta v_x = \frac{\delta\,\cos\gamma}{\cos(\alpha/2)}, \qquad
\Delta v_y = \frac{\delta\,\sin\gamma}{\cos(\alpha/2)}
\qquad\Rightarrow\qquad
\frac{\Delta v_y}{\Delta v_x} = \tan\gamma
$$

The split between $v_x$ and $v_y$ depends **only on $\gamma$** — $\alpha$ drops out of the
ratio entirely. $\alpha$ instead sets the *overall amplitude* via $1/\cos(\alpha/2)$, which
diverges as $\alpha\to180°$ (near-antiparallel headings — the typical same-platform
ascending/descending case) and is modest ($1/\cos(45°)=\sqrt2$) at $\alpha=90°$ (orthogonal
crossing). Near-antiparallel geometry amplifies a common-mode bias the most, regardless of
pixel position; pixel position only determines *which* of $v_x$/$v_y$ absorbs more of it.

### 2. Independent measurement noise (left panel: random $dr$ error)

Uncorrelated noise in $p_A$ and $p_D$ separately propagates via the **row norms** of
$\mathbf{A}$ — a different combination than case 1's row *sums*:

$$
\sigma_{v_x} = \sqrt{A_{00}^2+A_{01}^2}\,\sigma_p, \qquad
\sigma_{v_y} = \sqrt{A_{10}^2+A_{11}^2}\,\sigma_p
$$

At $\alpha=90°$ exactly, $\sin\alpha=1,\cos\alpha=0$ and $\mathbf{A}$ collapses to
$\begin{pmatrix}\cos\beta&-\sin\beta\\\sin\beta&\cos\beta\end{pmatrix}$ — an orthonormal
**rotation matrix**, so $\sigma_{v_x}=\sigma_{v_y}$ (for $\sigma_{p_A}=\sigma_{p_D}$) at
*every* pixel position: perfectly isotropic error propagation. Departing from $\alpha=90°$
introduces $\gamma$-dependent anisotropy — e.g. at $\alpha=170°$ (near-antiparallel), the
$\sigma_{v_x}/\sigma_{v_y}$ ratio ranges from $\approx0.12$ to $\approx8$ depending on
$\gamma$, vs. exactly $1$ at every $\gamma$ for $\alpha=90°$.

The left panel's y-axis converts a raw, per-measurement range-offset error $\sigma_{dr}$
(metres, over a `REPEAT_DAYS`-day repeat cycle) into $\sigma_p$ (m/yr) via
$\sigma_p = \frac{365.25}{\Delta t\,\sin\psi}\,\sigma_{dr}$ (the same scaling `mosaic3d` itself
applies to range offsets — see `mosaic3d.md` Step 2), so the plotted values are directly in
m/yr per metre of raw measurement noise, using each region's real incidence angle.

### Practical implication

Crossing geometry near $\alpha=90°$ (orthogonal headings) is doubly favorable — it minimizes
amplification of common-mode bias *and* gives uniform, isotropic sensitivity to independent
measurement noise, regardless of pixel position. Same-platform ascending/descending pairs
($\alpha$ near $180°$) are the worst case on both counts, though whether $v_x$ or $v_y$ is hit
harder by a given common-mode bias depends entirely on $\gamma$, not on $\alpha$ itself.

---

## Code

- [`nisarerrors/plotVerticalSensitivity.py`](../nisarerrors/plotVerticalSensitivity.py) — this
  script.
- `headingAndIncidence()` ports `common/computeHeading.c`'s lookDir-aware cross-track-heading
  formula onto real state vectors via `utilities.geodatrxa.geodatrxa.llzPtToRA()` (a faithful
  port of `llToImageNew.c`'s zero-Doppler Newton solve).
- `computeA()` is a direct translation of `common/initRoutines.c`'s `computeA()`.

## Further reading

- GrIMP's `mosaicSource/Documents/mosaic3d.md` — full `mosaic3d` algorithm documentation,
  including "Error Analysis: Crossing-Geometry Sensitivity" (Step 1) — the source derivation
  this page summarizes — and "Vertical-Motion Correction Sign Convention" (the $\cot\psi$
  derivation for common-mode vertical-motion errors).
- `mosaicSource/common/initRoutines.c` (`computeA`, `computeVxy`) and
  `mosaicSource/common/computeHeading.c` — the C source this script's geometry is verified
  against.
