# Tiled phase unwrapping of noisy NISAR interferograms: the cost of residues and a masking remedy

**Status: draft, 2026-11-14.** Based on real-data diagnostics for 5 Antarctica40 frames and a synthetic test of
the first 10 of 100 frames. The synthetic test is running on all 100 frames, and this note will be extended when
it finishes (see [Plan](#plan)).

## Summary

Unwrapping with snaphu split into 4 x 4 tiles, followed by a single-tile re-optimization (`-S`), was about 2.5
times faster than untiled snaphu on 80 MHz Greenland interferograms, with identical results. On 40 MHz
Antarctic interferograms the same approach was often slower than untiled snaphu. We show that the slowdown comes
from residue-dense noise, mostly open ocean and sea ice. The tiles solve the noise independently, they disagree
along their seams, and the single-threaded re-optimization then has to reconcile hundreds of thousands of
residues.

Removing the noise before unwrapping restores the benefit of tiling, but only if the noisy pixels are both masked
and filled with a smooth phase. Masking alone leaves the residues in the problem, and a crude fill creates new
ones. We use a mask built from the residues themselves: 15 x 15 pixel blocks whose residue density exceeds 0.25
per loop, restricted to regions at least three blocks thick.

In a synthetic test with a known answer, masked tiled unwrapping was 7-19 times faster than the unmasked untiled
unwrapping used in production (median 10 times) on the first 10 frames. It was as accurate on 8 of the 9 frames
that could be scored. On the remaining frame the mask, not the tiling, increased the error, though only slightly
in absolute terms (2.5% against 1.2% of pixels with the wrong cycle count), by disconnecting pieces of good data
whose relative cycle offsets then went wrong.

## 1. Introduction

snaphu solves phase unwrapping as a minimum-cost network flow problem. Its cost grows with the number of
**residues**: 2 x 2 pixel loops whose wrapped phase differences do not sum to zero, each of which must be
connected to a residue of the opposite sign or to the edge of the data. In coherent data residues are rare. In
decorrelated data, such as open water, about one loop in three is a residue.

To unwrap large interferograms faster, snaphu can divide the image into overlapping tiles, unwrap them in
parallel, and assemble the result. With `-S` it then re-optimizes the assembled solution as a single tile, which
removes the errors tiling would otherwise introduce at the seams. That last step is single-threaded and works on
the whole scene, but it is normally quick. It starts from the assembled tile solution, which on clean data is
already almost right, so the solver has little left to change.

We are evaluating this tiled mode (with a build of snaphu 2.0.7 that removes a fixed per-tile delay and speeds up
part of the network solver) for unwrapping NISAR RUNW products. The goal is to make unwrapping cheap enough to
redo routinely, for example after correcting the interferogram, without changing the answer. On 68 Greenland
frames the tiled mode was 2.5 times faster and gave identical results. On the first Antarctic coastal frame it
was 1.5 times *slower*. This note explains why, describes a preprocessing step that removes the penalty, and
tests whether that step costs accuracy.

## 2. Data and methods

### 2.1 Data

We use 100 Antarctica40 RUNW products: 40 MHz (`_4000_`), HH, unwrapped on the 13 x 16 look RUNW grid. They were
selected across tracks, pass directions, latitudes and acquisition dates, and classed as **coastal** (60 frames,
footprint touching ocean, ice shelf or the grounding zone) or **interior** (40 frames, fully grounded). Each
product's own unwrapped phase, rewrapped, together with its coherence, forms the interferogram we unwrap. JPL's
unwrapped phase and connected components serve as a reference solution.

### 2.2 Unwrapping

All runs use snaphu 2.0.7 with the statistical cost mode `SMOOTH`, an MCF initialization, `NCORRLOOKS 139`, and the
RUNW's 13 x 16 looks. The two modes compared are:

- **untiled**: the unmodified snaphu release, as used for production unwrapping;
- **tiled with `-S`**: the modified build, 4 x 4 tiles with 200-pixel overlaps, 16 processes, followed by the
  single-tile re-optimization. Its output is byte-identical to the unmodified release for the same mode.

Runs were made on a 48-core machine shared with other work. Wall-clock times are therefore indicative to perhaps
20%.

### 2.3 Comparing solutions

Two unwrapped solutions can differ by an arbitrary whole number of cycles between connected components. We
therefore compare them per *joint* component, meaning each pair of overlapping components from the two
solutions. Within each pair we remove the pair's median cycle offset and count the pixels that still differ by
one or more cycles. Against the real data this measures **disagreement** between two minimum-cost solutions, not
error.

### 2.4 Synthetic test

To measure error, we build a synthetic interferogram with a known answer for each frame:

- **Truth.** The phase predicted from the GrIMP velocity model on the RUNW grid (`simPhase/velSim.tif`).
- **Noise.** At each pixel, the phase error of an L-look interferogram with that pixel's coherence, drawn from the
  exact distribution of a two-channel complex Wishart matrix (two random draws per pixel). The noise therefore
  follows the frame's real coherence pattern.
- **Coherence.** The RUNW's own coherence, corrected in two ways:
  - **Low-coherence bias.** Estimated coherence reads high when the true coherence is low (open ocean reads about
    0.08 when the truth is near 0), so it is debiased for 139 looks.
  - **Fringe-slope bias.** It reads low where the phase changes quickly across the estimation window, because
    the phasors partly cancel. We divide by the expected loss, sinc(g_r/2) x sinc(g_a/2), computed from the
    truth's own phase gradients g_r and g_a in radians per pixel. Cells steeper than about 2.4 radians per pixel
    are left decorrelated, as they would be in reality.
- **Number of looks.** L is fitted per frame (between 16 and 208) so that the synthetic interferogram has the same
  residue density as the real one in each coherence range.

Each frame is unwrapped twice: **unmasked untiled**, the production approach, and **masked tiled `-S`**. For some frames, earlier diagnostic runs also provide the two mixed combinations, which
separate the effect of the mask from the effect of tiling. Each solution is scored as the share of pixels, in
its connected components, whose cycle count differs from the component's median offset against the truth. Both
solutions are scored over the same pixels: those the mask keeps.

## 3. Results

### 3.1 Tiling is slow when residue-dense noise remains in the interferogram

These runs use the real RUNW phase, rewrapped, with no mask. How much time tiling saves depends on how many
residues are left after the tiled pass, and on where they are.
Coherent frames have few residues and tile well. The interior frame `t085D_f129` (31,000 residues) was 3.5 times
faster tiled, in line with the Greenland results. Coastal frames carry large areas of open ocean, where nearly
every third loop is a residue. The coastal frame `t107D_f123` has 714,000 residues (8.7% of all loops). Untiled
snaphu took 775 s and the tiled mode 1,133 s, with identical results.

The tiled run averaged only 1.4 busy cores. Its parallel stage is short, and almost all of the time goes to the
single-threaded re-optimization. Figure 1 shows why that step is expensive. 77% of the residues lie in open
ocean, and they are concentrated in a few tiles (five of the sixteen hold 107,000-139,000 residues each, most of
the rest fewer than 5,000). Several seams run straight through the densest areas. Each tile pairs the residues
near its edges differently, and the re-optimization then has to reconcile them across the whole scene.

![Figure 1](images/snaphuTiling_t107D_residuesAndTiles.png)
*Figure 1. Residue density on `t107D_f123` (purple, log scale), the 4 x 4 snaphu tiles (orange; shading = 200-pixel
overlap) and the ocean more than 10 km from land or ice (blue). Numbers give the residues in each tile.*

JPL's production processing does not remove this noise either. The products list a water mask among their
inputs, but the recorded unwrapping preprocessing names only a sub-swath mask. The unwrapped phase over the
ocean is noise (median gradient 2.1 radians per pixel against 0.34 over ice), and coherent ocean pixels appear in
connected components. Unmasked untiled unwrapping is therefore the right reference for production.

### 3.2 Noise must be masked and filled

Removing the noise from snaphu's problem takes more than a mask. snaphu's byte mask makes masked pixels free to
cut through, but their wrapped phase still defines the residues to be solved:

- **Mask only.** Masking the ocean more than 10 km offshore in the byte mask alone left the untiled time unchanged
  (827 s against 775 s). The tiled time fell to 164 s, because the masked regions no longer contributed cost to
  the re-optimization.
- **Mask and fill.** Replacing the masked phase with a fill removes the masked residues altogether.

The kind of fill matters. A nearest-neighbour fill leaves a phase step where fills from two sides of a thin masked
band meet. The end of such a step acts as a residue, and on `t085D_f129` it produced a straight one-cycle cut
through coherent ice. We instead fill with a distance-weighted complex average of the surrounding phase (a
normalized Gaussian convolution of the unit phasors, widening until it finds data). That is close to the
`distance_interpolator` JPL uses, and it leaves almost no residues inside filled areas.

A second preprocessing detail concerns the edges of each frame. Where a sample is invalid in either radar image,
the RUNW stores a phase of exactly 0, about 5-10% of each frame along the swath edges. Treating that 0 as data
creates a phase step along the whole edge. We fill these pixels in the same way and keep them in the unwrapping,
as JPL does, but leave them out of all comparisons.

### 3.3 A residue-based mask

Rather than relying on a coastline, we derive the mask from the residues themselves. The residue density in
15 x 15 pixel blocks is close to two-valued: blocks are either nearly free of residues or close to the one in three
expected for random phase. We mask blocks above a threshold, but only in regions at least three blocks thick, so
that narrow low-coherence streaks across the ice are left to snaphu's coherence weighting. The mask is computed
from the data alone and needs neither a coastline nor geolocation.

Table 1 compares two thresholds on five frames. On the two interior frames the mask hardly matters. On the
coastal frames a threshold of 0.10 removes more residues and runs faster, while 0.25 removes only near-random
blocks and agrees slightly better with JPL. In every case the tiled and untiled solutions were identical, with
one exception. On `t138A_f141` at 0.25, 0.007% of pixels differed: with 415,000 residues left, the two methods
started their non-convex optimization from different points and reached slightly different minima.

*Table 1. Real-data runs with the residue mask at thresholds 0.10 and 0.25. Colour marks a 0.25 value that differs
substantially from the 0.10 value (time more than 1.5 times; disagreement changed by more than half, or by more
than 5 percentage points): red worse, green better.*

| Frame | Class | Masked (0.10 / 0.25) | Untiled (0.10 / 0.25) | Tiled + `-S` (0.10 / 0.25) | Disagreement with JPL, % (0.10 / 0.25) |
|---|---|---|---|---|---|
| `t048A_f137` | interior | 5.4% / 4.4% | 601 / 693 s | 82 / 104 s | 0.031 / 0.031 |
| `t085D_f129` | interior | 0.2% / 0.1% | 709 / 694 s | 238 / 217 s | 1.283 / 1.282 |
| `t107D_f123` | coastal | 27.5% / 23.7% | 746 / 855 s | 116 / 127 s | 0.329 / $\color{green}{\text{0.079}}$ |
| `t055D_f125` | coastal | 54.8% / 46.5% | 1,187 / 1,583 s | 99 / $\color{red}{\text{199}}$ s | 0.166 / 0.143 |
| `t138A_f141` | coastal | 57.5% / 39.3% | 1,571 / $\color{red}{\text{2,933}}$ s | 374 / $\color{red}{\text{1,427}}$ s | 24.1 / $\color{red}{\text{29.9}}$ |

Frames that are mostly noise are the difficult case (Figure 2). On `t138A_f141`, about 40-60% of the frame is
masked. The good data that remains is a set of fragments (glacier tongues, ice-shelf pieces, rough ground) linked
only by narrow strips. JPL's solution connected these fragments through rough terrain that our mask removes, and
our solution and JPL's disagree on the relative cycle offsets of several of them (24-30% of compared pixels).

Neither solution can be confirmed from the data. A local comparison of each disagreeing fragment against the
velocity-model phase favoured our solution on some fragments and JPL's on others, and was inconclusive where the
model itself is several cycles off, near fronts and glacier tongues. Raising the threshold kept more of the
connecting terrain but made both modes much slower, without reducing the disagreement.

![Figure 2](images/snaphuTiling_t138A_maskThresholds.png)
*Figure 2. The residue mask on `t138A_f141` at four thresholds (purple). Red marks masked pixels that JPL's solution
placed in a connected component.*

A carrier correction made no difference to unwrapping cost. Adding the GrIMP rubber-sheet carrier correction, with
or without the ionospheric correction, changed the residue count of `t085D_f129` by 1.6% and the tiled time by
less than the run-to-run spread. The correction is smooth, so it can only flip residues where the phase is
already noise.

### 3.4 Synthetic test: speed and accuracy

Table 2 gives the results for the first 10 frames. Masked tiled `-S` unwrapping is **7-19 times faster** than
unmasked untiled unwrapping, median 10.1 times. Together the ten frames took 11,200 s unmasked untiled and 923 s
masked tiled `-S`. The largest gains are on large, noisy coastal frames, where untiled unwrapping is
slowest.

*Table 2. Synthetic test, first 10 frames: unmasked untiled against masked tiled `-S`. "Wrong" is the share of
scored pixels whose cycle count is wrong. Red marks a masked tiled `-S` error substantially worse than the unmasked
untiled one (more than double, or more than 0.5 percentage points higher); no value is substantially better.*

| Frame | Class | Size (Mpx) | L | Residues | Masked | Unmasked untiled | Masked tiled `-S` | Speedup | Wrong, unmasked untiled | Wrong, masked tiled `-S` |
|---|---|---|---|---|---|---|---|---|---|---|
| `t003D_f121` | coastal | 7.0 | 16 | 1,915k | 77.8% | 770 s | 53 s | 14.6 | (no components) | (no components) |
| `t014A_f143` | coastal | 6.6 | 128 | 767k | 31.6% | 677 s | 81 s | 8.4 | 0.530% | 0.554% |
| `t048A_f137` | interior | 12.4 | 208 | 287k | 2.7% | 514 s | 72 s | 7.1 | 0.000% | 0.000% |
| `t055D_f125` | coastal | 9.3 | 64 | 1,638k | 47.4% | 1,543 s | 100 s | 15.4 | 0.061% | 0.068% |
| `t064A_f140` | coastal | 14.6 | 16 | 3,408k | 55.8% | 3,233 s | 168 s | 19.3 | 0.203% | 0.203% |
| `t078D_f122` | coastal | 3.9 | 16 | 989k | 73.9% | 270 s | 23 s | 11.6 | 0.228% | 0.227% |
| `t085D_f129` | interior | 12.5 | 208 | 54k | 0.1% | 769 s | 104 s | 7.4 | 0.000% | 0.000% |
| `t107D_f123` | coastal | 8.2 | 64 | 686k | 23.3% | 870 s | 102 s | 8.6 | 0.960% | 0.921% |
| `t144A_f143` | coastal | 10.5 | 64 | 1,420k | 26.1% | 2,332 s | 191 s | 12.2 | 0.970% | 0.929% |
| `t170A_f145` | coastal | 3.5 | 32 | 820k | 59.7% | 222 s | 30 s | 7.3 | 1.152% | $\color{red}{\text{2.488\%}}$ |

On 8 of the 9 frames that could be scored, the masked tiled `-S` error is within 0.03 percentage points of the
unmasked untiled error (median 0.215% against 0.216%). One frame, `t003D_f121`, is almost entirely noise, and neither run
produced a connected component. The exception is `t170A_f145`, where the masked tiled `-S` error is twice the unmasked untiled error.

Diagnostic runs on that frame separate the two changes. Unmasked tiled `-S` unwrapping gives exactly the
unmasked untiled error (1.152%), and masked untiled unwrapping gives exactly the masked tiled `-S` error (2.488%). The extra error therefore comes
from the mask, as on `t138A_f141` in the real data, and not from tiling. On every frame where both were run, tiled
and untiled unwrapping of the same input gave the same error.

Most of the remaining error is not caused by noise. On `t107D_f123` nearly all of it lies in a narrow rim along the
ice front (Figure 3). There the velocity-model phase jumps by up to 38 radians per pixel between ice and ocean,
which no unwrapper can recover from a wrapped interferogram, and the masked and unmasked runs make the same
mistakes. In real data such fronts lose coherence before the jump matters, so this part of the synthetic error is
pessimistic.

![Figure 3](images/snaphuTiling_t107D_syntheticErrors.png)
*Figure 3. Synthetic test on `t107D_f123`. Top: phase gradient of the truth; blue outlines where it exceeds π
radians per pixel. Middle and bottom: pixels with the wrong cycle count in the masked and unmasked solutions
(red too high, blue too low; grey masked). From the first synthetic run, before the fringe-slope correction; the
corrected run's error rates differ by less than 0.04 percentage points.*

## 4. Discussion

Tiling with `-S` never changed the unwrapped answer in a way that mattered. Its speed, however, depends on the
residues left after the tiled pass. snaphu's tiled pass is parallel, but the re-optimization that makes the result
seam-free is not. On clean data the assembled tiles are almost right and that step converges quickly. When the
tiles meet in dense noise, each tile pairs the noise residues differently, the assembled solution is far from the
whole-scene optimum, and the re-optimization can take longer than unwrapping the whole scene untiled. Removing the
noise first leaves the assembled solution nearly right again, so almost all of the work stays in the parallel stage. Untiled unwrapping time
depends mostly on frame size and changes little when the noise is masked (Table 1), so the gain from masking goes
almost entirely to the tiled mode.

The mask carries a risk that tiling does not. On frames where good data survives only as fragments, masking the
noise between them removes the paths that tied their phases together. Their relative cycle offsets are then fixed
only through masked or narrow regions, and can come out wrong. This is the mechanism behind the extra error on
`t170A_f145` and the disagreement with JPL on `t138A_f141`. It is a property of the mask, and a gentler rule on land
(for example, masking only residue regions connected to open ocean) may avoid it. In production, the per-component
ambiguity correction against the velocity model, which can treat separated pieces of a component independently,
would also limit the damage.

Agreement with JPL's solution is a weak test. Both it and ours are minimum-cost solutions to the same wrapped
phase, and their disagreements are confined to fragments whose offsets the data do not determine. The synthetic
test is the more informative measure, with two limitations. Its truth is a velocity model, which is aliased at
some ice fronts. And its noise is calibrated to each frame's residue statistics, not to its actual phase.

## 5. Conclusions

1. Tiled unwrapping with `-S` gives the same result as untiled unwrapping. On noisy coastal frames it is fast only
   if the residue-dense noise is removed first.
2. Removing the noise requires both masking and a smooth fill of the masked pixels, together with filling the
   zero-valued sub-swath edges of the RUNW.
3. A residue-density mask (15 x 15 blocks, threshold 0.25, regions at least three blocks thick) makes masked tiled
   unwrapping 7-19 times faster than the production approach on the first 10 synthetic frames, with the same
   accuracy on 8 of 9 scorable frames.
4. Where the mask does cost accuracy, the cost so far is small (2.5% against 1.2% of pixels with the wrong cycle
   count, on one frame), and it comes from disconnecting fragments of good data. The mask rule, not the use of
   tiles, is what needs refining for such frames.

## Plan

When all 100 synthetic frames are done:

1. **Full results.** Extend Table 2 and summarize speed and accuracy by class (coastal or interior), size and
   noise level.
2. **Frames that lose accuracy.** Check each one as `t170A_f145` was checked (masked untiled and unmasked tiled),
   to attribute the loss to the mask or to tiling.
3. **Gentler masks.** Test mask variants on those frames, for example masking only residue regions connected to
   open ocean, or a different thickness rule.
4. **Production.** Settle the mask rule, and add the preprocessing (residue mask, smooth fill of masked and
   sub-swath pixels) ahead of tiled unwrapping in production.
5. **Real-data check.** Repeat the real-data comparison on a sample of frames with the final preprocessing.

## Appendix: files

Work directory `/Volumes/insar4/ian/Data/NISAR/Antarctica/unwrapTest40/`:

| Path | Contents |
|---|---|
| `frames.csv` | the 100 frames and their class |
| `bin/runUnwrapTest.py` | the original unmasked comparison (stopped after 2 frames) |
| `diagOceanMask/diagResidueMask2.py` | real-data runs with the residue mask; results in `residueMask2/`, `residueMask2_t0.25/` and `residueMask2carrier*/` |
| `diagOceanMask/patchVsSynthetic.py` | local test of disagreeing fragments against the velocity-model phase |
| `diagOceanMask/synthTest.py` | the synthetic test; results in `synth2/{none,mask0.25}/<frame>/result.json` |
| `diagOceanMask/synth2/queueAll.sh` | the queue working through all 100 frames |
