# Tiled phase unwrapping of noisy NISAR interferograms: the cost of residues and a masking remedy

**Status: draft, 2026-11-15.** Based on real-data diagnostics for 5 Antarctica40 frames and a synthetic test. The
first version of the synthetic test (Sections 3.4 to 3.7, Appendix A) used a truth that is zero, not missing, where
the velocity model has no coverage. A corrected version is running on 200 frames; the Summary gives its provisional
results, and the results sections will be redone from it (see [Plan](#plan)).

## Summary

These experiments are designed to see if the wall-clock time for snaphu can be reduced by using tiled rather than
untiled unwrapping. Specifically, we investigate tiling with the `-S` option, using a slightly modified version of
snaphu 2.0.7 (Section 2.2). Tiling splits the interferogram into tiles (here 4 x 4) that are unwrapped in parallel
and then assembled; with `-S`, snaphu then re-optimizes the assembled solution as a single tile so that the result
has no errors along the tile seams. Production unwrapping of NISAR RUNW products uses
untiled snaphu, which runs on a single core because its network solver is sequential.

On 80 MHz Greenland interferograms, tiled unwrapping with `-S` was about 2.5 times faster than untiled unwrapping,
with identical results. On 40 MHz Antarctic interferograms it was often slower. The slowdown comes from
residue-dense noise, mostly open ocean and sea ice. The tiles solve the noise independently, they disagree along
their seams, and the single-threaded re-optimization then has to reconcile hundreds of thousands of residues.

Removing the noise before unwrapping restores the benefit of tiling, but only if the noisy pixels are both masked
and filled with a smooth phase. Masking alone leaves the residues in the problem, and a crude fill creates new
ones. We use a mask built from the residues themselves: 15 x 15 pixel blocks whose residue density exceeds 0.25
per loop, restricted to regions at least three blocks thick.

In a synthetic test with a known answer, masked tiled unwrapping is a median 8.1 times faster than the unmasked
untiled unwrapping used in production (provisional, 119 of 200 frames, 50 of them coastal). Scaled to the 1,803
frame locations of the Antarctica40 project (30% coastal), one pass of unwrapping would take about 80 h instead of
about 570 h, about 7 times less. Accuracy is unchanged on almost every frame: on 103 of the 116 frames that could
be scored the two error rates differ by less than 0.07 percentage points, and the mean share of pixels with the
wrong cycle count is 1.2% in both. The mask makes a few frames substantially more accurate (5) and a few less
accurate (3); the losses come from disconnecting pieces of good data whose relative cycle offsets then go wrong.

On about one frame in fifteen, snaphu's assembly of the tiles fails to converge. The `-S` re-optimization repairs
the result, so the answer is unaffected, but these frames lose most of the speedup and take about half of the
masked tiled time. Wider tile overlaps reduce this on the worst frames but cost about as much on the others.

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

We performed experiments both with actual phase data and with a larger set of synthetic interferograms based on a
reference velocity map. The two answer different questions. The actual data show how snaphu behaves on real
noise, including open ocean, sea ice and decorrelated shear margins, and how long it takes. But they have no
known answer, so two solutions can only be compared with each other or with JPL's (Section 2.3). The synthetic
interferograms combine the phase predicted by the velocity map with noise matched to each frame's real coherence
(Section 2.4). They have a known answer, so they measure error as well as time, on many more frames.

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
- **tiled with `-S`**: a slightly modified build of snaphu 2.0.7 (below), 4 x 4 tiles with 200-pixel overlaps, 16
  processes, followed by the single-tile re-optimization.

The modified build changes 125 of snaphu's roughly 20,200 lines of C (114 lines added, 11 removed, in four
files). None of the changes alter the unwrapped result: wherever the unmodified release runs to completion, the
output is byte-identical. We checked the first two changes against the unmodified release on CoSSC data and on 68
Greenland frames, and the third against the build without it on Antarctic frames here.
The changes are:

- **No delay between tile processes** (1 line). snaphu pauses before starting each tile process, in proportion to
  the tile size. The pause is unnecessary, because each tile writes its own files and the parent waits for them.
- **Faster ancestor search** (94 lines added, 6 removed). The tree solver repeatedly searches for the common
  ancestor of two nodes (`FindApex`) by walking up the tree one node at a time. The modified build caches a pointer
  to an ancestor every 32 levels and invalidates the cache when the tree changes, so the search skips most of the
  walk.
- **No abort on large flows next to masked pixels** (19 lines added, 4 removed). snaphu aborts if any arc carries
  more flow than its 16-bit cost arithmetic allows: about 160 cycles at the default `NSHORTCYCLE` of 200. A tiled solution can leave that much flow on
  arcs next to masked pixels, which carry no cost and which the solver already ignores when it tracks the largest
  flow. The modified build applies the check to the other arcs only (Section 3.6).

Runs were made on machines shared with other work: a 48-core machine and, for the corrected synthetic test, also a
96-core machine. Both runs of a frame were always made on the same machine. Wall-clock times are therefore
indicative to perhaps 20%, and speedups are more reliable than absolute times.

### 2.3 Comparing solutions

Two unwrapped solutions can differ by an arbitrary whole number of cycles between connected components. We
therefore compare them per *joint* component, meaning each pair of overlapping components from the two
solutions. Within each pair we remove the pair's median cycle offset and count the pixels that still differ by
one or more cycles. Against the real data this measures **disagreement** between two minimum-cost solutions, not
error. In the synthetic test the truth takes the place of the second solution, and the same count measures error
(Section 2.4).

Specifically, we evaluate intra-component error, not inter-component error. A component whose phase is off from the
reference by a whole number of cycles counts as correct; only pixels whose cycle count differs from the rest of
their component count as wrong. Offsets between components are left to later processing, such as the
per-component ambiguity correction against the velocity model.

### 2.4 Synthetic test

To measure error, we build a synthetic interferogram with a known answer for each frame:

- **Truth.** The phase predicted from the GrIMP velocity model on the RUNW grid (`simPhase/velSim.tif`).
- **Noise.** At each pixel, the phase error of an L-look interferogram with that pixel's coherence, drawn from the
  exact distribution of a two-channel complex Wishart matrix (two random draws per pixel). The noise therefore
  follows the frame's real coherence pattern.
- **Coherence.** The RUNW's own coherence, corrected in two ways:
  - **Low-coherence bias.** Coherence is estimated as the magnitude of an average of N products of the two images'
    complex samples, normalized by their powers. When the true coherence is zero, those products have random
    phases, and the magnitude of their average is not zero but about 1/sqrt(N), because random phasors do not cancel
    exactly. For the N = 139 effective looks of the RUNW coherence (13 x 16 samples, fewer independent ones), this
    floor is about 0.085, and open ocean, where the true coherence is near zero, reads about 0.08. The bias shrinks
    as the true coherence rises. We remove it with the approximation E[g_est^2] = g^2 + (1 - g^2)/N, i.e.
    g = sqrt((N g_est^2 - 1)/(N - 1)), set to zero below the floor. Without this correction, every decorrelated
    pixel would get a little spurious coherence, and the synthetic noise would be too weak exactly where it matters
    most.
  - **Fringe-slope bias.** It reads low where the phase changes quickly across the estimation window, because
    the phasors partly cancel. We divide by the expected loss, sinc(g_r/2) x sinc(g_a/2), computed from the
    truth's own phase gradients g_r and g_a in radians per pixel. Cells steeper than about 2.4 radians per pixel
    are left decorrelated, as they would be in reality.
- **Number of looks.** The RUNW is formed from 13 x 16 = 208 samples, but we do not use 208 as L. Neighbouring
  samples are correlated, so the number of independent looks is smaller, and the real phase noise also departs from
  the ideal model in ways the coherence does not capture: sea ice and shelf fronts change between acquisitions, and
  the two coherence corrections above are approximate. With L = 208 the synthetic interferogram would be cleaner
  than the real one on some frames and noisier on others. What makes a frame hard to unwrap is its residues, so we
  choose L for each frame, from 16, 32, 64, 128 and 208, to make the synthetic interferogram's residue density in
  each coherence range (bins from 0 to 1) match the real interferogram's as closely as possible. L is therefore a
  calibration constant that sets how hard the synthetic frame is to unwrap, not a physical look count. The fitted
  values span the whole range: coastal frames often need a low L (their real phase is noisier than their coherence
  implies), and interior frames usually fit best at 208.

### 2.5 Experiment protocol

Each synthetic frame is unwrapped twice: **unmasked untiled**, the production approach, and **masked tiled `-S`**.
Both runs use the same noisy interferogram (the noise is seeded per frame), and the masked run additionally applies
the residue mask and smooth fill (Section 3.3). Each solution is scored as the share of pixels, in its connected
components, whose cycle count differs from the component's median offset against the truth (Section 2.3). Both
solutions are scored over the same pixels: those the mask keeps.

Where the masked tiled run is substantially less accurate, we also unwrap the masked input untiled with the
unmodified snaphu. Comparing that run with the other two separates the effect of the mask from the effect of
tiling.

Pixels where the truth is exactly zero are treated as missing: the velocity model stores zero, not a no-data
value, where it has no coverage, which is mostly over ice shelves, ice tongues and the ocean. Those pixels are
filled and left out of the unwrapping and the scoring, like pixels with zero coherence or in sub-swath gaps. A
frame is skipped if the model covers less than half of its coherent pixels (RUNW coherence above 0.5), and is
replaced by another frame of the same class chosen by the same rule (Section 2.1). The first version of the test
treated zero truth as data (see the status note at the top).

Each machine runs at most five untiled and three tiled unwrappings at a time, and both runs of a frame are made on
the same machine. Wall-clock and CPU times come from `/usr/bin/time`.

## 3. Results

We first use the actual phase data to show why tiling alone does not speed up unwrapping on noisy Antarctic frames
(Section 3.1), and to develop a preprocessing step that fixes it: masking the residue-dense noise and filling it
smoothly (Sections 3.2 and 3.3). We then use the synthetic interferograms to measure what that approach gains in
time and costs in accuracy (Sections 3.4 and 3.5). Finally, we describe two problems in snaphu's tiled mode that
the tests exposed, an overflow abort and a slow assembly of the tiles on some frames (Sections 3.6 and 3.7).

### 3.1 Why tiling with no masking does not work

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
overlap) and the ocean more than 10 km from land or ice (blue), the mask of our first attempt to remove the noise
(Section 3.2). Numbers give the residues in each tile.*

JPL's production processing does not remove this noise either. The products list a water mask among their
inputs, but the recorded unwrapping preprocessing names only a sub-swath mask. The unwrapped phase over the
ocean is noise (median gradient 2.1 radians per pixel against 0.34 over ice), and coherent ocean pixels appear in
connected components. Unmasked untiled unwrapping is therefore the right reference for production.

### 3.2 Noise must be masked and filled

Our first attempt to remove the noise was an ocean mask, the approach JPL's processing appears to anticipate (it
lists a water mask among its inputs). We masked the ocean more than 10 km from land or ice, taken from our
Antarctic ice mask and mapped into radar geometry with the frame's orbit (blue in Figure 1, 27% of `t107D_f123`).
On that frame it worked for the tiled mode, but it is not the approach we kept:

- It removes only open ocean. Sea ice, the 10 km coastal margin, fjords and decorrelated areas on land keep their
  noise, and on many coastal frames these hold much of it.
- It depends on an external mask, which must be current and geocoded into radar geometry for every frame.
- The residue mask of Section 3.3, which needs only the interferogram itself, did better on the same frame (116 to
  127 s tiled, against 164 s).

The ocean-mask runs also showed that removing the noise from snaphu's problem takes more than a mask. snaphu's
byte mask makes masked pixels free to cut through, but their wrapped phase still defines the residues to be solved:

- **Mask only.** Masking the ocean in the byte mask alone left the untiled time unchanged (827 s against 775 s).
  The tiled time fell from 1,133 s to 164 s, because the masked regions no longer contributed cost to the
  re-optimization.
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

### 3.5 Synthetic test: all 100 frames

Table 3 and Figure 4 summarize the 100 frames (60 coastal, 40 interior); Appendix A lists every frame. Masked tiled
`-S` unwrapping took a median of 97 s per frame against 670 s for unmasked untiled unwrapping (Figure 4a). The median
speedup is 7.3 times (Figure 4b), and the total time fell 5.9 times, from 52.6 h to 9.0 h.

**Accuracy.** Most frames are unchanged (Figure 4c, d). On 69 of the 94 frames that could be scored the two errors
are within 0.07 percentage points of each other, and on 81 within 0.5 points. The remaining differences go both
ways: ten frames are substantially better and six substantially worse (Table 3). The gains are large on a few
frames (`t151A_f138` 12.8% to 5.9%, `t044D_f123` 3.6% to 0.09%, `t097D_f125` 1.8% to 0.6%), and the losses are
large on one (`t032D_f121` 6.5% to 32.9%). Over all scored frames the mean error rose from 3.6% to 3.8% and the
median fell from 0.098% to 0.074%; the losses (36 percentage points summed) outweigh the gains (17 points), most of
them on `t032D_f121`.

All six losses come from the mask, not from tiling: unwrapping the masked input untiled gives the same error to
within 0.02 percentage points on every one of them (`t170A_f145`, `t059D_f122`, `t032D_f121`, `t162A_f142`,
`t123A_f138` and `t111A_f144`). The mask works in both directions. It helps where it removes noise that untiled
snaphu would otherwise route through good data, and it hurts where it cuts good data into pieces whose relative
cycle offsets the remaining data no longer fix (Section 4).

*Table 3. Synthetic test, all 100 frames, by class. Totals are summed wall-clock times. "Wrong" is the median share
of scored pixels with the wrong cycle count. "Substantially worse" and "substantially better" use the rule of
Table 2 in each direction.*

| Set | Frames | Unmasked untiled, total | Masked tiled `-S`, total | Speedup, median (range) | Wrong, unmasked untiled | Wrong, masked tiled `-S` | Substantially worse | Substantially better |
|---|---|---|---|---|---|---|---|---|
| all | 100 | 52.6 h | 9.0 h | 7.3 (0.8-141.8) | 0.098% | 0.074% | 6 of 94 | 10 of 94 |
| coastal | 60 | 35.5 h | 5.6 h | 8.1 (0.8-41.2) | 0.965% | 0.854% | 6 of 54 | 7 of 54 |
| interior | 40 | 17.1 h | 3.4 h | 6.6 (0.8-141.8) | 0.000% | 0.000% | 0 of 40 | 3 of 40 |

![Figure 4](images/snaphuTiling_syntheticSummary.png)
*Figure 4. Synthetic test, all 100 frames. (a) Wall time of each run, with medians (dashed) and totals; runs longer
than 4,000 s are in the last bin (unmasked untiled | masked tiled `-S`). (b) Speedup, stacked by class, with the
median (dashed); speedups above 30 are in the last bin, with their values. (c) Share of scored pixels with the wrong
cycle count, masked tiled `-S` against unmasked untiled, on log axes (zero plotted at 0.001%); frames substantially
worse are labelled. (d) Change in that share, masked tiled `-S` minus unmasked untiled; changes beyond ±3 points
are in the end bins, with their values.*

### 3.6 An overflow abort

Three frames (`t059D_f122`, `t116D_f126`, `t081D_f127`) first failed. snaphu aborted at the start of the `-S`
re-optimization with "((Maximum flow) * NSHORTCYCLE) too large". Each tile leaves its masked pixels at an arbitrary
cycle count, so the assembled solution jumped by up to 245 cycles across masked areas, while snaphu allows at most
160 (32,000 divided by `NSHORTCYCLE` 200). Arcs next to masked pixels carry no cost, and snaphu's solver loop
already ignores them when it tracks the largest flow, so the modified build now applies the initial check to the
other arcs only. On frames that passed the check its output is byte-identical to the previous build. One frame of
the second set (`t060A_f139`) still fails, because a jump of 389 cycles lies between unmasked pixels (Section 3.7).

### 3.7 Slow frames: the assembly of the tiles

snaphu assembles the tiles by solving a secondary network for the cycle offset of each piece of each tile. On 14 of
the 100 frames this solver stopped early ("No overall cost reduction for too many iterations"), and the assembled
solution contains wrong offsets of tens of cycles between pieces (127 cycles across a 180-pixel boundary on
`t081D_f127`), which the `-S` pass has to repair. The answer is repaired, but the time is lost: these 14 frames
have a median speedup of 3.0 against 7.8 for the other 86, and take two thirds of the masked tiled `-S` time
(6.1 of 9.0 h). If they ran at the median speedup of the others, the total would fall to about 5.6 h, a speedup of
about 9.4 rather than 5.9.

The cause lies in the tiles, not in the solver. The solver works with costs clipped to 16-bit integers, but an exact
solution of the same secondary network for `t081D_f127` with 64-bit costs is only 16% cheaper and still contains a
96-cycle offset. Masked and decorrelated areas cut tiles into pieces whose unwrapped solutions disagree in the tile
overlaps, and the secondary costs then favour large offsets. Changes to the assembly did not help: scaling the
secondary costs differently, solving the secondary network twice at two cost resolutions, and choosing the offsets
between whole tiles jointly from all overlaps (which raised the agreement in the overlaps from 50% to 85% of
pixels).

Wider tile overlaps give the same answer at about the same total cost. We reran masked tiled `-S` with overlaps of
400 pixels instead of 200 on 17 frames: 13 slow frames (the ones whose assembly stopped early, and three that were
slow without stopping) and 4 normal ones. Every frame's error was unchanged to within 0.004 percentage points. The
slow frames mostly ran faster (`t081D_f127` 3,193 s to 1,466 s, `t173D_f127` 866 s to 428 s, `t138A_f141` 1,736 s to
1,108 s), and the assembly converged on six of them, but `t111A_f144` ran slower (2,375 s to 4,494 s). Three of the
four normal frames ran about 20% slower, because the tiles are larger. In total the 17 frames took 12% less time
(12,516 s against 14,224 s), but they are mostly slow frames; weighting the two groups as in the 100-frame set
(Section 3.5) gives roughly the same total time either way. The one frame that aborted at 200 pixels
(`t060A_f139`, Section 3.6) runs at 400 pixels. These runs used the uncorrected synthetic inputs (Plan), which does
not affect the comparison of the two overlaps on the same input. With 3 x 3 tiles, `t081D_f127` took 1,924 s.

## 4. Discussion

Tiling with `-S` never changed the unwrapped answer in a way that mattered. Its speed, however, depends on the
residues left after the tiled pass. snaphu's tiled pass is parallel, but the re-optimization that makes the result
seam-free is not. On clean data the assembled tiles are almost right and that step converges quickly. When the
tiles meet in dense noise, each tile pairs the noise residues differently, the assembled solution is far from the
whole-scene optimum, and the re-optimization can take longer than unwrapping the whole scene untiled. Removing the
noise first leaves the assembled solution nearly right again, so almost all of the work stays in the parallel stage. Untiled unwrapping time
depends mostly on frame size and changes little when the noise is masked (Table 1), so the gain from masking goes
almost entirely to the tiled mode.

The mask cuts both ways. On frames where good data survives only as fragments, masking the noise between them
removes the paths that tied their phases together. Their relative cycle offsets are then fixed only through masked
or narrow regions, and can come out wrong. This is the mechanism behind the six frames the mask made worse, and
the disagreement with JPL on `t138A_f141`. On other frames the same mask removes noise that untiled snaphu would
otherwise connect through the data, and the result improves. A gentler rule on land (for example, masking only
residue regions connected to open ocean) may keep the gains and avoid the losses. In production, the per-component
ambiguity correction against the velocity model, which can treat separated pieces of a component independently,
would also limit the damage.

The fragmentation that hurts accuracy also costs time. When masked and decorrelated areas cut a tile into pieces,
the tiles disagree in their overlaps, snaphu's assembly fails to converge, and the `-S` pass does the work instead
(Section 3.7). Wider overlaps reduce this on the worst frames, but cost about as much on the others.

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
   unwrapping a median 7.3 times faster than the production approach, and 5.9 times faster in total, on 100
   synthetic frames, with the same accuracy on most of them.
4. The mask changes accuracy on about one frame in six, more often for the better (10 frames) than for the worse
   (6), but the losses are larger, and on balance it costs a little (a mean of 3.8% against 3.6% of pixels with the
   wrong cycle count). The losses come from disconnecting fragments of good data. The mask rule, not the use of
   tiles, is what needs refining.
5. On one frame in seven, snaphu's assembly of the tiles fails to converge and `-S` repairs it at a large cost in
   time. This is not a solver limitation. Wider tile overlaps (400 pixels) shift time from these frames to the
   others without changing the answer, and save little overall.

## Plan

Done for the first 100 frames: the full results (Section 3.5, Appendix A), and the attribution of every frame that
lost accuracy (all to the mask). Next:

1. **Second 100 frames.** Repeat the summary for the second set (P05023 products only), running now, and for both
   sets together.
2. **Tile overlaps.** Repeat the 400-pixel overlap comparison (Section 3.7) on a random sample of the corrected
   frames, together with a timing of stock snaphu against the modified build.
3. **Gentler masks.** Compare frames the mask helped (`t151A_f138`, `t044D_f123`) with frames it hurt (`t032D_f121`,
   `t172D_f127`), and test mask variants, for example masking only residue regions connected to open ocean, or a
   different thickness rule.
4. **Production.** Settle the mask rule and the tiling parameters, and add the preprocessing (residue mask, smooth
   fill of masked and sub-swath pixels) ahead of tiled unwrapping in production.
5. **Real-data check.** Repeat the real-data comparison on a sample of frames with the final preprocessing.

## Appendix A: synthetic test, all 100 frames

*Table A1. Synthetic test, the first 100 frames. Columns as in Table 2. Red marks a masked tiled `-S` error
substantially worse than the unmasked untiled one, green one substantially better (more than double or half, or
more than 0.5 percentage points apart).*

| Frame | Class | Size (Mpx) | L | Residues | Masked | Unmasked untiled | Masked tiled `-S` | Speedup | Wrong, unmasked untiled | Wrong, masked tiled `-S` |
|---|---|---|---|---|---|---|---|---|---|---|
| `t001A_f138` | interior | 14.0 | 208 | 8k | 0.0% | 485 s | 55 s | 8.9 | 0.000% | 0.000% |
| `t003D_f121` | coastal | 7.0 | 16 | 1,915k | 77.8% | 770 s | 53 s | 14.6 | (no components) | (no components) |
| `t004D_f133` | interior | 9.6 | 208 | 4k | 0.0% | 347 s | 28 s | 12.5 | 0.000% | 0.000% |
| `t007D_f123` | coastal | 8.6 | 64 | 1,202k | 40.7% | 499 s | 85 s | 5.9 | 0.771% | 0.788% |
| `t008D_f133` | interior | 6.9 | 208 | 553k | 1.4% | 303 s | 38 s | 7.9 | 0.000% | 0.000% |
| `t011D_f126` | coastal | 14.3 | 64 | 2,386k | 34.9% | 7,411 s | 180 s | 41.2 | 0.346% | 0.254% |
| `t011D_f128` | interior | 13.9 | 208 | 46k | 0.0% | 360 s | 54 s | 6.6 | 0.004% | 0.004% |
| `t013A_f143` | interior | 3.8 | 208 | 8k | 0.0% | 59 s | 11 s | 5.1 | 0.000% | 0.000% |
| `t013A_f144` | coastal | 14.7 | 64 | 2,563k | 41.8% | 4,053 s | 186 s | 21.8 | 0.911% | 0.924% |
| `t014A_f143` | coastal | 6.6 | 128 | 767k | 31.6% | 677 s | 81 s | 8.4 | 0.530% | 0.554% |
| `t017D_f123` | interior | 9.3 | 208 | 1,873k | 36.5% | 507 s | 82 s | 6.2 | 0.002% | 0.002% |
| `t018D_f131` | interior | 14.0 | 208 | 94k | 0.0% | 336 s | 67 s | 5.0 | 0.000% | 0.000% |
| `t020A_f139` | interior | 14.7 | 208 | 66k | 0.0% | 439 s | 67 s | 6.6 | 0.000% | 0.000% |
| `t024A_f143` | interior | 6.9 | 208 | 296k | 1.4% | 235 s | 41 s | 5.7 | 0.000% | 0.000% |
| `t024D_f119` | coastal | 4.7 | 16 | 1,387k | 88.1% | 304 s | 29 s | 10.5 | (no components) | (no components) |
| `t025D_f126` | interior | 6.2 | 208 | 534k | 8.2% | 308 s | 37 s | 8.3 | 0.000% | 0.000% |
| `t027A_f141` | interior | 14.8 | 128 | 152k | 0.0% | 487 s | 196 s | 2.5 | 0.001% | 0.001% |
| `t030D_f122` | coastal | 4.2 | 64 | 800k | 45.0% | 343 s | 51 s | 6.7 | 1.973% | $\color{green}{\text{1.459\%}}$ |
| `t032D_f121` | coastal | 13.6 | 32 | 2,355k | 45.5% | 1,403 s | 544 s | 2.6 | 6.512% | $\color{red}{\text{32.918\%}}$ |
| `t036A_f138` | interior | 14.4 | 128 | 25k | 0.0% | 830 s | 76 s | 10.9 | 0.000% | 0.000% |
| `t036A_f139` | coastal | 16.4 | 32 | 2,423k | 43.4% | 3,360 s | 290 s | 11.6 | 3.403% | 3.616% |
| `t037A_f142` | coastal | 10.0 | 32 | 1,608k | 37.0% | 664 s | 103 s | 6.4 | 0.117% | 0.130% |
| `t038A_f144` | coastal | 16.1 | 128 | 782k | 7.5% | 1,020 s | 136 s | 7.5 | 0.178% | 0.194% |
| `t044A_f139` | interior | 14.7 | 64 | 0k | 0.0% | 311 s | 39 s | 7.9 | 0.000% | 0.000% |
| `t044D_f123` | coastal | 11.2 | 64 | 2,195k | 51.3% | 1,365 s | 118 s | 11.6 | 3.553% | $\color{green}{\text{0.087\%}}$ |
| `t048A_f137` | interior | 12.4 | 208 | 287k | 2.7% | 514 s | 72 s | 7.1 | 0.000% | 0.000% |
| `t050A_f140` | coastal | 11.6 | 32 | 2,904k | 71.3% | 1,002 s | 115 s | 8.7 | 0.665% | 0.670% |
| `t052A_f144` | coastal | 5.0 | 16 | 1,053k | 60.0% | 738 s | 121 s | 6.1 | 17.528% | 17.463% |
| `t053D_f119` | coastal | 3.2 | 16 | 762k | 70.7% | 156 s | 20 s | 7.8 | (no components) | (no components) |
| `t055D_f125` | coastal | 9.3 | 64 | 1,638k | 47.4% | 1,543 s | 100 s | 15.4 | 0.061% | 0.068% |
| `t056A_f145` | coastal | 7.7 | 32 | 1,346k | 48.1% | 443 s | 73 s | 6.1 | 0.058% | 0.059% |
| `t057D_f127` | coastal | 7.3 | 64 | 695k | 28.0% | 689 s | 781 s | 0.9 | 26.289% | $\color{green}{\text{25.474\%}}$ |
| `t059D_f122` | coastal | 12.7 | 64 | 2,344k | 52.0% | 2,046 s | 162 s | 12.6 | 1.669% | $\color{red}{\text{3.669\%}}$ |
| `t059D_f123` | interior | 13.9 | 208 | 65k | 0.1% | 601 s | 89 s | 6.7 | 0.015% | $\color{green}{\text{0.000\%}}$ |
| `t064A_f140` | coastal | 14.6 | 16 | 3,408k | 55.8% | 3,233 s | 168 s | 19.3 | 0.203% | 0.203% |
| `t066D_f127` | coastal | 10.8 | 16 | 796k | 20.6% | 1,387 s | 216 s | 6.4 | 7.053% | 6.730% |
| `t067A_f145` | coastal | 5.1 | 32 | 420k | 18.1% | 442 s | 40 s | 11.0 | 0.167% | 0.303% |
| `t069D_f125` | coastal | 15.1 | 16 | 3,579k | 69.6% | 1,718 s | 138 s | 12.5 | (no components) | (no components) |
| `t071A_f136` | interior | 13.3 | 208 | 0k | 0.0% | 355 s | 62 s | 5.7 | 0.000% | 0.000% |
| `t072D_f127` | coastal | 13.2 | 208 | 907k | 18.9% | 1,062 s | 292 s | 3.6 | 8.730% | 8.757% |
| `t073D_f123` | coastal | 15.8 | 128 | 1,983k | 24.8% | 1,373 s | 246 s | 5.6 | 0.150% | $\color{green}{\text{0.020\%}}$ |
| `t074D_f121` | coastal | 8.2 | 32 | 1,852k | 56.5% | 890 s | 71 s | 12.5 | 0.069% | 0.064% |
| `t075A_f146` | coastal | 15.1 | 16 | 4,520k | 90.5% | 1,405 s | 127 s | 11.1 | (no components) | (no components) |
| `t075D_f123` | interior | 3.8 | 208 | 28k | 0.0% | 132 s | 17 s | 7.8 | 0.000% | 0.000% |
| `t076A_f142` | coastal | 2.4 | 16 | 725k | 92.1% | 105 s | 16 s | 6.5 | (no components) | (no components) |
| `t078D_f122` | coastal | 3.9 | 16 | 989k | 73.9% | 270 s | 23 s | 11.6 | 0.228% | 0.227% |
| `t080A_f139` | coastal | 10.1 | 208 | 599k | 12.1% | 582 s | 153 s | 3.8 | 13.538% | 13.553% |
| `t080D_f133` | interior | 10.0 | 208 | 408k | 2.8% | 328 s | 65 s | 5.0 | 0.000% | 0.000% |
| `t081D_f127` | coastal | 15.0 | 32 | 1,443k | 27.2% | 11,098 s | 3,193 s | 3.5 | 4.230% | 4.183% |
| `t082D_f130` | interior | 4.2 | 208 | 0k | 0.0% | 125 s | 13 s | 9.4 | 0.000% | 0.000% |
| `t083D_f125` | coastal | 14.7 | 32 | 1,840k | 27.3% | 1,232 s | 171 s | 7.2 | 0.041% | 0.041% |
| `t085D_f126` | coastal | 10.9 | 64 | 1,862k | 39.1% | 1,951 s | 114 s | 17.1 | 0.029% | 0.038% |
| `t085D_f129` | interior | 12.5 | 208 | 54k | 0.1% | 769 s | 104 s | 7.4 | 0.000% | 0.000% |
| `t086A_f136` | interior | 14.0 | 208 | 6k | 0.0% | 396 s | 49 s | 8.1 | 0.000% | 0.000% |
| `t087A_f142` | interior | 14.4 | 64 | 483k | 5.1% | 20,277 s | 143 s | 141.8 | 0.036% | $\color{green}{\text{0.016\%}}$ |
| `t087D_f125` | coastal | 14.7 | 32 | 1,287k | 17.6% | 1,704 s | 596 s | 2.9 | 3.473% | 3.470% |
| `t095A_f143` | coastal | 7.8 | 128 | 1,099k | 32.6% | 1,252 s | 84 s | 14.9 | 0.188% | 0.193% |
| `t097D_f125` | coastal | 15.8 | 128 | 4,036k | 70.7% | 5,958 s | 320 s | 18.6 | 1.779% | $\color{green}{\text{0.562\%}}$ |
| `t099D_f126` | coastal | 14.4 | 32 | 3,259k | 58.0% | 3,885 s | 144 s | 26.9 | 0.063% | 0.063% |
| `t102D_f124` | interior | 14.0 | 208 | 377k | 0.1% | 456 s | 89 s | 5.1 | 0.000% | 0.000% |
| `t103D_f127` | interior | 9.7 | 208 | 68k | 0.0% | 284 s | 48 s | 5.9 | 0.000% | 0.000% |
| `t106A_f140` | coastal | 14.3 | 208 | 2,217k | 22.3% | 5,412 s | 216 s | 25.1 | 1.607% | 1.165% |
| `t107D_f123` | coastal | 8.2 | 64 | 686k | 23.3% | 870 s | 102 s | 8.6 | 0.960% | 0.921% |
| `t108A_f139` | coastal | 15.4 | 32 | 2,277k | 40.4% | 15,157 s | 1,351 s | 11.2 | 49.517% | 49.389% |
| `t108D_f127` | interior | 14.7 | 208 | 1k | 0.0% | 474 s | 41 s | 11.5 | 0.000% | 0.000% |
| `t111A_f144` | coastal | 14.7 | 128 | 1,681k | 18.8% | 10,872 s | 2,375 s | 4.6 | 5.262% | $\color{red}{\text{7.684\%}}$ |
| `t113D_f126` | interior | 14.3 | 128 | 3,252k | 60.6% | 2,999 s | 175 s | 17.1 | 0.001% | 0.001% |
| `t115D_f127` | coastal | 10.8 | 208 | 854k | 20.0% | 590 s | 143 s | 4.1 | 13.704% | $\color{green}{\text{12.162\%}}$ |
| `t115D_f129` | interior | 12.0 | 208 | 14k | 0.0% | 507 s | 41 s | 12.4 | 0.027% | 0.027% |
| `t116A_f141` | interior | 5.4 | 208 | 7k | 0.0% | 130 s | 21 s | 6.1 | 0.000% | 0.000% |
| `t116D_f126` | coastal | 14.3 | 128 | 547k | 5.3% | 2,356 s | 261 s | 9.0 | 0.059% | 0.067% |
| `t119A_f140` | interior | 7.3 | 128 | 31k | 0.0% | 429 s | 33 s | 13.0 | 0.001% | 0.001% |
| `t121D_f124` | interior | 3.8 | 208 | 134k | 0.0% | 88 s | 16 s | 5.7 | 0.001% | 0.001% |
| `t123A_f138` | coastal | 10.8 | 208 | 60k | 0.1% | 622 s | 74 s | 8.4 | 3.055% | $\color{red}{\text{4.980\%}}$ |
| `t124A_f139` | coastal | 15.0 | 128 | 1,244k | 11.9% | 1,000 s | 165 s | 6.1 | 0.014% | 0.014% |
| `t126A_f140` | interior | 14.4 | 208 | 74k | 0.0% | 540 s | 95 s | 5.7 | 0.000% | 0.000% |
| `t127A_f142` | interior | 3.8 | 208 | 145k | 0.0% | 106 s | 16 s | 6.6 | 0.001% | 0.001% |
| `t128A_f145` | coastal | 7.8 | 32 | 1,274k | 47.5% | 642 s | 64 s | 10.1 | 0.183% | 0.219% |
| `t131D_f122` | coastal | 6.6 | 32 | 1,709k | 74.2% | 689 s | 55 s | 12.6 | 10.248% | $\color{green}{\text{8.281\%}}$ |
| `t136A_f140` | coastal | 14.3 | 32 | 3,656k | 69.4% | 4,146 s | 168 s | 24.7 | 0.079% | 0.079% |
| `t138A_f141` | coastal | 9.3 | 32 | 1,653k | 39.3% | 3,195 s | 1,736 s | 1.8 | 16.090% | 16.434% |
| `t144A_f135` | interior | 12.8 | 208 | 5k | 0.0% | 603 s | 49 s | 12.3 | 0.000% | 0.000% |
| `t144A_f143` | coastal | 10.5 | 64 | 1,420k | 26.1% | 2,332 s | 191 s | 12.2 | 0.970% | 0.929% |
| `t145D_f123` | coastal | 14.4 | 208 | 2,288k | 35.8% | 4,442 s | 213 s | 20.8 | 0.606% | 0.548% |
| `t146D_f123` | interior | 6.6 | 208 | 490k | 2.5% | 225 s | 36 s | 6.3 | 0.000% | 0.000% |
| `t148A_f142` | coastal | 5.4 | 128 | 1,413k | 77.2% | 319 s | 50 s | 6.4 | 0.283% | 0.274% |
| `t151A_f138` | interior | 13.9 | 32 | 157k | 0.7% | 20,095 s | 4,818 s | 4.2 | 12.760% | $\color{green}{\text{5.886\%}}$ |
| `t152A_f139` | coastal | 5.4 | 208 | 201k | 8.3% | 246 s | 297 s | 0.8 | 17.008% | 17.008% |
| `t162A_f142` | coastal | 14.4 | 16 | 938k | 11.2% | 1,510 s | 842 s | 1.8 | 14.030% | $\color{red}{\text{16.098\%}}$ |
| `t165D_f132` | interior | 14.1 | 208 | 0k | 0.0% | 541 s | 40 s | 13.5 | 0.000% | 0.000% |
| `t166A_f135` | interior | 12.5 | 128 | 120k | 0.0% | 3,541 s | 4,504 s | 0.8 | 14.547% | 14.547% |
| `t166A_f139` | coastal | 7.0 | 64 | 848k | 35.2% | 1,326 s | 556 s | 2.4 | 49.623% | 49.287% |
| `t168A_f144` | coastal | 14.3 | 208 | 440k | 3.1% | 595 s | 111 s | 5.4 | 0.001% | 0.002% |
| `t169D_f131` | interior | 3.8 | 208 | 299k | 8.3% | 154 s | 86 s | 1.8 | 3.914% | 3.913% |
| `t170A_f145` | coastal | 3.5 | 32 | 820k | 59.7% | 222 s | 30 s | 7.3 | 1.152% | $\color{red}{\text{2.488\%}}$ |
| `t171A_f141` | interior | 15.1 | 128 | 118k | 0.2% | 1,620 s | 674 s | 2.4 | 0.976% | 0.625% |
| `t171A_f144` | coastal | 11.9 | 32 | 799k | 11.2% | 804 s | 110 s | 7.3 | 0.030% | 0.033% |
| `t172D_f127` | coastal | 5.1 | 208 | 426k | 19.4% | 1,607 s | 712 s | 2.3 | 10.280% | 10.504% |
| `t172D_f132` | interior | 6.6 | 208 | 75k | 0.0% | 213 s | 40 s | 5.3 | 0.001% | 0.001% |
| `t173D_f127` | coastal | 12.7 | 128 | 887k | 16.1% | 965 s | 866 s | 1.1 | 8.000% | 8.184% |

## Appendix B: files

Work directory `/Volumes/insar4/ian/Data/NISAR/Antarctica/unwrapTest40/`:

| Path | Contents |
|---|---|
| `frames.csv` | the 100 frames and their class |
| `bin/runUnwrapTest.py` | the original unmasked comparison (stopped after 2 frames) |
| `diagOceanMask/diagResidueMask2.py` | real-data runs with the residue mask; results in `residueMask2/`, `residueMask2_t0.25/` and `residueMask2carrier*/` |
| `diagOceanMask/patchVsSynthetic.py` | local test of disagreeing fragments against the velocity-model phase |
| `diagOceanMask/synthTest.py` | the synthetic test; results in `synth2/{none,mask0.25}/<frame>/result.json` |
| `diagOceanMask/synth2/queueAll.sh` | the queue working through all 100 frames |
| `diagOceanMask/synthSummary.py` | summary table, per-frame table and figure for the synthetic test |
| `bin/snaphu-2.0.7-grimpopt-maskflow`, `bin/snaphu-2.0.7-grimp-maskflow.patch` | the modified snaphu build used for masked tiled `-S` since 2026-11-14, and its patch against snaphu 2.0.7 |
| `diagOceanMask/reanchorTest/attribute.sh` | masked untiled run of one frame, to separate the effect of the mask from that of tiling |
| `diagOceanMask/reanchorTest/ovl400.sh`, `scoreOvl400.py` | the 400-pixel overlap test (`reanchorTest/ovl400/`), kept apart from the main runs |
| `diagOceanMask/synth2/queue2.sh`, `frames2.csv`, `bin/selectMoreFrames.py` | the second 100 frames |
