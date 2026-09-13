# validateMosaic — the error-calibration harness

Acceptance test for a velocity mosaic's **formal errors**. Built to calibrate the crossing-pair
over-count parameter `rho` (`mosaicSource/Documents/mosaic3d.md`, §"Interpreting the formal
errors"), and generalised so any solver can be checked the same way — in particular a replacement
crossing-orbit method such as the per-pixel normal-equation design in
`mosaicSource/Documents/crossingOrbitRedundancy.md` §5.

## What it measures

| symbol | meaning |
|---|---|
| `e` | per-pixel formal error, the `ex`/`ey` bands |
| `d` | `S1 − NISAR` difference at a pixel, `vx` and `vy` pooled |
| `RMS(e)` | `sqrt(mean(e²))` — the population prediction from the formal errors |
| `std(d)` | observed standard deviation of the differences |
| **`k`** | **`std(d)/RMS(e)`. 1.0 is calibrated; `k > 1` means reported errors are too small** |

**Use `RMS(e)`, never `mean(e)` or `median(e)`.** Variances average, standard deviations do not.
An early version of this analysis used a median and biased `k` high by 25–40%; a second version
compared a MAD-based scale against `RMS(e)` and produced a spurious "calibrated at rho = 0" result.
Both are recorded in the report below so they are not repeated.

## Where everything lives

| piece | path |
|---|---|
| S1 reference cache (25 MB, with YAML provenance sidecar) | `/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject/Notebooks/calValReport/data/Greenland/S1reference/` |
| test-mosaic root, one directory per run | `/Volumes/insar1/.../Release/velocity/errCal/` |
| merged products the tool reads | `.../errCal/merged/` |
| NISAR input list (symlinked as `track-all` in each run dir) | `/Volumes/insar1/.../Release/masterInput/` |
| full experimental record, including discarded runs | `.../errCal/report/errCalResults.md` |

The reference is the mean of 9 NSIDC-0725 annual mosaics (2014-12-01 → 2023-11-30) at
`overviewLevel=2`, built by `makeS1Reference.py` (see `makeS1Reference.md`). Its inter-annual sigma
indicates ~1 m/yr or better over stable ground, so a several-m/yr difference is NISAR's, not S1's.

## Running it

```bash
python -m nisarerrors.validateMosaic opt_t10000_phase opt_t37_offsets opt_t10000_both --coverage
```

Tags are merged-product prefixes resolved against `--mergedDir`, or full paths. A prefix `P` needs
`P_vx.tif`, `P_vy.tif`, `P_ex.tif`, `P_ey.tif`.

Statistics are taken on the footprint **common to every product in the call**, so results are
comparable across runs — but they therefore shift slightly depending on which products are compared
together. Always validate a new solver against a baseline **in the same invocation**.

## Baseline to beat

The shipped configuration (`rho = 0.6`), at the production thresholds
(`timePhaseThresh 10000`, `timeThresh 37`), 409,161 px:

| product | RMS(e) | std(d) | k |
|---|---|---|---|
| phase | 2.985 | 2.985 | **1.00** |
| offsets | 10.578 | 8.854 | 0.84 |
| both | 3.257 | 3.588 | 1.10 |

Coverage at the same setting — the scale is right but the shape is not Gaussian:

| product | <1σ | <2σ | <3σ | 90% | 95% | 99% | 99.9% |
|---|---|---|---|---|---|---|---|
| phase | 86.4% | 96.3% | 98.3% | 1.19 | 1.71 | 3.89 | 9.66 |
| both | 83.6% | 95.9% | 98.3% | 1.30 | 1.82 | 3.75 | 8.46 |
| *Gaussian* | *68.3%* | *95.4%* | *99.7%* | *1.65* | *1.96* | *2.58* | *3.29* |

## Adding a new test run

Copy an existing run directory (they are self-contained) and edit only the flags under test:

```bash
cd /Volumes/insar1/.../Release/velocity/errCal
cp -r opt_t10000_phase myNewMethod          # keeps template, mosaicsSetup.yaml, track-all symlink
# edit myNewMethod/mosaic.template.yaml: regionID (must be unique, hyphens only) and baseFlags
cd myNewMethod && makemosaic.py --nThreads 10 --interval multiYear \
    --firstdate 2025-06-01 --lastdate 2026-11-30 \
    --mosaicsSetupFile mosaicsSetup.yaml --template mosaic.template.yaml --noLabel
```

Then merge the 48 tiles into the `merged/` prefix the tool reads:

```bash
T=myNewMethod/Vel-2025-06-01.2026-11-30/tiff
for pair in "vx vx" "vy vy" "vv v" "ex ex" "ey ey"; do
  set -- $pair
  gdalbuildvrt -q -overwrite merged/myNewMethod.$1.vrt $T/mosaic-*.$2.tif
  gdal_translate -q -of GTiff -co COMPRESS=DEFLATE -co TILED=YES \
      merged/myNewMethod.$1.vrt merged/myNewMethod_$1.tif
done
```

### Grid, and a trap that invalidates everything

1600 m grid, `948 × 1712`, an exact 8× decimation of the 200 m release grid sharing its NW corner,
so it aligns with the S1 reference at `overviewLevel=2`:

```
sx = 948   sy = 1712   dx = dy = 1600   xll = -658300.0   yll = -3377500.0   nx = 6   ny = 8
```

**`-fl` is in output-grid PIXELS, not metres.** Production `-fl 10` at 200 m is a 2 km feather;
carried unchanged to a 1600 m grid it becomes 16 km and changes inter-frame blending — and hence
the error field under test — far more than anything you are trying to measure. **All runs on this
grid must use `-fl 1`.** Every template here already does; preserve it when copying.

## Notes for validating a new crossing-orbit solver

- **The threshold experiment does not transfer.** `rho` was calibrated by raising the crossing time
  threshold, which adds *pairs* while the contributing *images* stay fixed (image count moves 0.5%
  from T=12 to T=37 while pair count moves 2.5×), so the measured error is flat by construction and
  a correctly scaled formal error must be flat too. A pairless method has no such knob. The
  analogous lever is the number of contributing **images** (vary the date range), where the error
  genuinely *should* fall as `1/sqrt(n)`; the test is whether the reported error tracks that.
- **Disable the over-count correction on any pairless path.** `inflatePairOverCount()` would
  double-correct a result that has no over-counting by construction. Run with `-noPairOverCount`,
  or gate the call.
- **Expect `Cov = N⁻¹` to be optimistic.** It assumes independent measurements. The closest
  experimental analogue run here — a sandbox build consuming each image at most once per pixel,
  which has no combinatorial over-counting at all — still needed `rho ≈ 0.50–0.58` to calibrate.
  Extrapolating those runs to `rho = 0` (pure independence) gives `k ≈ 2.25` (phase) and `2.35`
  (offsets), i.e. **reported errors ~2.3× too small**. This is an extrapolation, not a measured
  point — the nearest actual run is `rho = 0.3` — but the linearity of `RMS(e)²` in `rho` is
  verified to 0.2%. The physical cause is that date-pairs sharing an acquisition inherit that
  frame's baseline residual and ionospheric screen, so `Cov(d)` is not diagonal. The principled fix
  is a low-rank shared-error term, `Cov = D + Σ_f σ_f² u_f u_fᵀ` with one `u_f` per frame, inverted
  via Woodbury and used as the GLS weight.
- **Outliers matter more than the tables suggest.** For crossing offsets in `vy` the worst 1% of
  pixels hold ~48% of the variance. A joint solve discards the accidental robustness the pair
  scheme gets from averaging many pairs, so plan for the IRLS pass rather than treating it as
  optional. Phase has the milder tail, which is why it is the right path to convert first.
