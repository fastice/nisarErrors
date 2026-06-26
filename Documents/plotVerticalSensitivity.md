# plotVerticalSensitivity — Reference

This page is a reference for how `mosaic3d`'s crossing-orbit velocity inversion responds to two
qualitatively different error sources — uncompensated vertical motion, and independent
range/phase measurement noise — as a function of crossing geometry and pixel position. It
documents the derivation, the figure, and the findings for both InSAR phase and speckle-tracked
range/azimuth offsets. It is not a how-to-run guide — see the bottom of the page for the CLI
invocation.

---

## 1. Setup

`mosaic3d`'s crossing-orbit inversion (`common/initRoutines.c`'s `computeA`/`computeVxy`,
documented in full in `mosaicSource/Documents/mosaic3d.md` §"Full Inversion") solves

$$
\begin{pmatrix} v_x \\ v_y \end{pmatrix} =
\left(\mathbf{I} - \mathbf{A}\mathbf{B}\right)^{-1} \mathbf{A}
\begin{pmatrix} p_A \\ p_D \end{pmatrix}, \qquad
\mathbf{D} = (\mathbf{I}-\mathbf{AB})^{-1}\mathbf{A}
$$

where $p_A$, $p_D$ are the ascending/descending LOS measurements scaled to horizontal-velocity
units (m/yr), and $\mathbf{B}$ is the DEM-slope correction. This page (and `plotVerticalSensitivity.py`)
uses the **flat-terrain approximation** ($\mathbf{B}=0$, i.e. $\mathbf{D}=\mathbf{A}$), since the
question here is the pure crossing-geometry sensitivity, not DEM-slope effects:

$$
\mathbf{A} = \frac{1}{\sin^2\!\alpha}
\begin{pmatrix}
\cos\beta - \cos\alpha\cos(\alpha+\beta) & \cos(\alpha+\beta) - \cos\alpha\cos\beta \\
\sin\beta - \cos\alpha\sin(\alpha+\beta) & \sin(\alpha+\beta) - \cos\alpha\sin\beta
\end{pmatrix}
$$

with $\alpha = H_A - H_D$ (heading difference) and $\beta = \phi - H_A$ ($\phi$ = pixel azimuth
in polar-stereographic coordinates, $H_A$ = ascending heading). For the sensitivity analysis
below it's more natural to define $H_{\text{mean}}=(H_A+H_D)/2$ and
$\gamma = \phi - H_{\text{mean}}$ (pixel azimuth relative to the *mean* heading, rather than to
$H_A$ alone) — full derivation in `mosaic3d.md` §"Error Analysis: Crossing-Geometry Sensitivity".

$H_A$, $H_D$ are the **cross-track** heading at the pixel (`common/computeHeading.c`), not the
along-track flight direction — these differ by a curved-Earth, incidence-angle-dependent amount
(confirmed ~6°–16° for this dataset), not a clean ±90°, so flight-direction heading cannot be
substituted. `headingAndIncidence()` ports `computeHeading.c`'s lookDir-aware formula onto real
orbital state vectors via `utilities.geodatrxa.geodatrxa.llzPtToRA()` (a verified line-for-line
port of `llToImageNew.c`'s zero-Doppler Newton solve).

Three real ascending/descending crossing pairs (northern/central/southern Greenland) provide the
example geometries, with $H_A$, $H_D$, $\alpha$, incidence angle $\psi$, and the real pixel
position $\gamma_{\text{real}}$ all computed at runtime from real `geodat26x16.geojson` files
(not assumed/idealized values):

| Region | $\alpha$ | $\psi_{\text{mean}}$ | $\gamma_{\text{real}}$ |
|---|---|---|---|
| North | 100.8° | 37.8° | 288.3° |
| Central | 106.4° | 37.7° | 271.2° |
| South | 123.4° | 37.9° | 268.6° |

---

## 2. The two error sources, and how each sensitivity is calculated

### 2a. Common-mode error (e.g. uncompensated vertical motion)

A common-mode error is the **same physical quantity contributing to both** $p_A$ and $p_D$. The
leading example is uncompensated vertical motion $v_z^{\text{SMB}}$ (unmodeled tide or
submergence/emergence correction — see `mosaic3d.md` §"Vertical-Motion Correction Sign
Convention"): a vertical rate $v_z^{\text{SMB}}$ projects onto each LOS via that image's own
incidence angle, contributing
$\delta_A = v_z^{\text{SMB}}\cot\psi_A$, $\delta_D = v_z^{\text{SMB}}\cot\psi_D$ to $p_A, p_D$
respectively (**not** simply $v_z^{\text{SMB}}\times1$ — that only holds at $\psi=45°$).

For equal contributions $\delta$ to both measurements (the case plotted: $\delta=1$ m/yr of
uncompensated $v_z^{\text{SMB}}$, scaled per-region by each image's real incidence angle via
$\cot\psi_A,\cot\psi_D$), $\mathbf{A}\cdot(\delta_A,\delta_D)^T$ reduces via sum-to-product
identities to a clean closed form:

$$
\Delta v_x = \frac{\delta\,\cos\gamma}{\cos(\alpha/2)}, \qquad
\Delta v_y = \frac{\delta\,\sin\gamma}{\cos(\alpha/2)}
\qquad\Rightarrow\qquad
\frac{\Delta v_y}{\Delta v_x} = \tan\gamma
$$

**Key result:** the split between $v_x$/$v_y$ depends *only* on $\gamma$ (pixel position relative
to mean heading) — $\alpha$ drops out of the ratio entirely. $\alpha$ instead sets the *overall
amplitude* via $1/\cos(\alpha/2)$, which diverges as $\alpha\to180°$ (near-antiparallel headings —
the typical same-platform ascending/descending case) and is modest ($1/\cos45°=\sqrt2$) at
$\alpha=90°$ (orthogonal crossing). Near-antiparallel geometry amplifies a common-mode bias the
most, regardless of where the pixel sits; pixel position only determines *which* of $v_x$/$v_y$
absorbs more of it.

In code: `sensitivities()` computes this directly as `A @ [cotPsiA, cotPsiD]` rather than via the
closed form, so it is exact (not a small-angle or other approximation) for any $\alpha,\beta$.

### 2b. Independent (uncorrelated) error in $p_A$, $p_D$ separately

This is the case of random, uncorrelated measurement noise in each LOS separately
($\sigma_{p_A},\sigma_{p_D}$ independent — e.g. random range/phase measurement noise). Variance
propagates as the diagonal of the output covariance:

$$
\sigma_{v_x}^2 = A_{00}^2\,\sigma_{p_A}^2 + A_{01}^2\,\sigma_{p_D}^2, \qquad
\sigma_{v_y}^2 = A_{10}^2\,\sigma_{p_A}^2 + A_{11}^2\,\sigma_{p_D}^2
$$

i.e. for $\sigma_{p_A}=\sigma_{p_D}$, $\sigma_{v_x},\sigma_{v_y}$ are the **row norms** of
$\mathbf{A}$ — a fundamentally different combination of $\mathbf{A}$'s entries than case (2a)'s
row *sums*, which is why the two cases have qualitatively different $\gamma$/$\alpha$ dependence.

**Key result:** at $\alpha=90°$ exactly, $\sin\alpha=1,\cos\alpha=0$ and $\mathbf{A}$ collapses to
$\begin{pmatrix}\cos\beta&-\sin\beta\\\sin\beta&\cos\beta\end{pmatrix}$ — an orthonormal
**rotation matrix** — so $\sigma_{v_x}=\sigma_{v_y}$ at *every* pixel position: perfectly
isotropic error propagation. Departing from $\alpha=90°$ introduces $\gamma$-dependent
anisotropy — e.g. at $\alpha=170°$ (near-antiparallel), the $\sigma_{v_x}/\sigma_{v_y}$ ratio
ranges from $\approx0.12$ to $\approx8$ depending on $\gamma$, vs. exactly 1 at every $\gamma$
for $\alpha=90°$.

The figure's left column converts a raw, per-measurement range-offset error $\sigma_{dr}$
(metres, over a `REPEAT_DAYS`-day repeat cycle) into $\sigma_p$ (m/yr) via
$\sigma_p = \frac{365.25}{\Delta t\,\sin\psi}\,\sigma_{dr}$ — the same scaling `mosaic3d` itself
applies to range offsets (`mosaic3d.md` Step 2) — using each region's real incidence angle.

**Practical implication of both results together:** crossing geometry near $\alpha=90°$ is
*doubly* favorable — it minimizes amplification of common-mode bias *and* gives uniform,
isotropic sensitivity to independent measurement noise, regardless of pixel position.
Same-platform ascending/descending pairs ($\alpha$ near 180°) are the worst case on both counts.
The three real Greenland crossing geometries above ($\alpha=101°$–$123°$) sit meaningfully closer
to the favorable 90° case than to the unfavorable 180° case.

---

## Figure — both sensitivities vs. pixel position, for three real crossing geometries

![vertical sensitivity](images/verticalSensitivity.png)

Left column: independent-error sensitivity $\sigma_{v_x},\sigma_{v_y}$ (m/yr, per metre of raw,
independent, equal range-offset/phase-equivalent noise in each LOS, for a 12-day repeat — see
§3 below for what "metre of raw error" means differently for phase vs. offsets). Right column:
common-mode response $\Delta v_x,\Delta v_y$ per 1 m/yr of uncompensated $v_z^{\text{SMB}}$.
Both swept over pixel position $\gamma=\phi-H_{\text{mean}}$ (0°–360°), with a secondary top axis
showing actual longitude and a dashed vertical line marking each region's real pixel position
$\gamma_{\text{real}}$ relative to its real crossing point.

At the real crossing points (all three regions land near $\gamma_{\text{real}}\approx270°$–288°),
the common-mode panel shows $\Delta v_y$ near its extremum and $\Delta v_x$ near zero — i.e. for
these specific real geometries, uncompensated vertical motion biases $v_y$ much more than $v_x$.
This was confirmed against a real `mosaic3d` run with an applied vertical correction, which
showed exactly this pattern ($v_x$ insensitive, $v_y$ sensitive) — resolving an earlier apparent
discrepancy that traced back to two now-fixed bugs in this script's example geometry (wrong
heading-computation method; an invalid ascending/descending track pairing for the North region),
not to any bug in `computeA`/`computeVxy` itself.

---

## 3. Findings for phase vs. range/azimuth offsets

**The geometric sensitivity formulas above (§2a, §2b) are identical for phase and offsets** —
confirmed directly in `mosaic3d.md`: Step 2 (`make3DOffsets`, range offsets) "[constructs] the
same $\mathbf{A}$ and $\mathbf{B}$ matrices (same as Step 1)... The crossing-geometry error
sensitivity... applies identically here — same $\mathbf{A}$, same $\gamma$/$\alpha$ dependence."
Unlike the squint analysis (`Documents/plotSquintError.md`), where phase has a genuinely separate
geometric exposure (the baseline-decomposition step), there is no such asymmetry here — both
measurement types feed the same $p_A, p_D \to (v_x,v_y)$ machinery via the identical
$\mathbf{A}$ matrix.

What *does* differ between phase and offsets is the **scaling from raw measurement to
$p_A$/$p_D$**, and — much more importantly in practice — the **typical magnitude of the raw
measurement noise itself**.

### Phase

$$
p_i = \frac{365.25}{\frac{4\pi}{\lambda_i}\,N_{\text{days},i}\,\sin\psi_i}\,\phi_i
$$

The input noise $\sigma_\phi$ combines baseline-parameter covariance (propagated through the 6×6
covariance matrix from the baseline fit) and tiepoint noise $\sigma_{\text{tp}}$
(`computePhiZM3d`, `mosaic3d.md` §"Phase Error"):

$$
\sigma_\phi = \sqrt{\mathbf{v}^T \mathbf{C}\, \mathbf{v} + \min(\pi,\,\sigma_{\text{tp}})^2}
$$

Converting $\sigma_\phi$ to an equivalent raw range-noise (the "metre of $dr$" unit this page's
figure is plotted in) is exactly $\sigma_{\delta r} = \frac{\lambda}{4\pi}\sigma_\phi$. At
L-band ($\lambda\approx0.24$ m), even a generous $\sigma_\phi\sim0.5$ rad of unwrapped-phase
noise corresponds to only $\sigma_{\delta r}\sim1$ cm of equivalent range noise — phase is an
intrinsically very precise measurement of $\delta r$. (This is an illustrative order-of-magnitude
figure, not a verified measurement from a specific real `mosaic3d` run.)

### Range/azimuth offsets (speckle tracking)

$$
p_A = \frac{365.25\,\Delta r_A}{\Delta t_A \sin\psi_A}
$$

The input noise combines the interpolated offset-field uncertainty, DEM-induced range error, and
baseline-parameter uncertainty (`mosaic3d.md` §"Step 2", point 6):

$$
\sigma_R = \sqrt{\sigma_{\text{off}}^2 + \sigma_{\text{dem}}^2 + \sigma_{\text{base}}^2}
$$

Speckle-tracked offset noise is set by cross-correlation precision (chip size, correlation,
oversampling — NISAR ATBD §5.4), which is typically dm-to-metre scale, not cm scale.

### Practical implication

Because the *geometric amplification factor* (row-norm of $\mathbf{A}$, plotted in the figure's
left column) is identical for both measurement types, but phase's typical equivalent range noise
is roughly one to two orders of magnitude smaller than offsets' typical range-offset noise, the
**same crossing-geometry sensitivity translates into a much smaller absolute velocity error for
phase than for offsets** at the same pixel/geometry. This is consistent with why `mosaic3d`'s
error-weighted combination (`mosaic3d.md` §"Error-Weighted Combination") generally favors phase
over offsets where both are available and coherence permits phase unwrapping — not because the
geometry treats them differently, but because phase is simply a quieter measurement of the same
underlying quantity.

---

## Usage

```
plotVerticalSensitivity [--projectDir DIR] [-o OUTPUT.png] [--show]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--projectDir DIR` | `/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject` | Root directory containing the example geodats. |
| `-o, --output PATH` | `verticalSensitivity.png` | Output plot path. |
| `--show` | off | Display interactively in addition to saving. |

The per-region track/frame relative paths (`REGIONS` in the script) are GIT64/NISAR-Greenland
example-specific and not currently overridable from the command line — edit the script directly
for a different project's example pairs.

## Code

[`nisarerrors/plotVerticalSensitivity.py`](../nisarerrors/plotVerticalSensitivity.py)
