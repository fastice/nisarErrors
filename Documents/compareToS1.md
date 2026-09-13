# compareToS1

Compares a NISAR velocity mosaic against a cached Sentinel-1 reference and reports whether
`mosaic3d`'s formal errors (`ex`/`ey`) are a correct description of the observed differences.

Build the reference first — see [makeS1Reference](makeS1Reference.md).

```bash
compareToS1 --overviewLevel 2 --tag stage0 \
  --product ".../all-phase-squint/.../GL-all-p-sq_..._*_v01.0" --label phase \
  --product ".../roffsets-squint/.../GL-all-r-sq_..._*_v01.0"  --label offsets \
  --product ".../all-squint/.../GL-all-rp-sq_..._*_v01.0"      --label both
```

## Why the obvious comparison is not enough

The natural thing to compute is `std(S1 - NISAR)` against `mean(ex)`. That was the starting
point of this investigation and it gave, over stable ground (S1 speed < 50 m/yr):

| product | S1 diff x/y | formal x/y |
|---|---|---|
| range offsets | 9.8 / 7.5 | 0.88 / 0.86 |
| phase | 3.0 / 3.6 | 1.02 / 1.04 |
| both | 4.0 / 7.3 | 0.76 / 0.85 |

Two things are wrong with reading those columns against each other.

**The statistics are not matched.** `std` is outlier-sensitive; `mean`/`median` of the formal
error is not. Every one of these difference distributions has `std/MAD` between 2.5 and 3.1
(a Gaussian gives 1.0), so the left column is inflated by a minority of large-difference pixels.
Compared robustly, the offsets are ~6x under-scaled, not ~11x.

**A single ratio cannot say what is wrong.** "The sigmas are globally too small" and "the sigmas
are fine except in bad areas" produce the same number and need different fixes.

So the tool reports robust and non-robust scales side by side, and adds the two statistics that
do discriminate.

## The statistics

### kGlobal

```
kGlobal = robustScale(difference) / median(formal)
```

`1.0` means calibrated. `6.0` means the formal errors are six times too small. Reported both
robustly (`1.4826*MAD`) and as `std`, so the effect of tails is visible rather than baked in.

All scales are also reported **after removing the median difference**. A one-signed bias is not
something a formal error bar models, so including it overstates the miscalibration.

### kDecile — the one that decides the fix

`k` recomputed within each decile of the formal error.

- **Flat** → the formal error is a good *relative* predictor and is simply scaled wrong. One
  global factor fixes it.
- **Declining** (large `k` where the formal error is smallest) → the signature of a **missing
  term added in quadrature**. A fixed extra error dominates a small sigma and barely moves a
  large one, so the ratio falls as the formal error grows. The fix is to add the missing term,
  not to rescale.
- **Rising** → the formal error is not even a good relative predictor; something per-pixel is
  missing.

Measured on the released products, this is what separated the two error models:

```
phase   vx  formal 0.56 0.63 0.71 0.78 0.86 0.94 1.04 1.15 1.35 2.04
            k      1.42 1.38 1.81 1.58 1.22 1.36 0.91 0.96 0.91 1.30   <- flat, ~1
offsets vx  formal 0.33 0.39 0.45 0.52 0.61 0.71 0.83 1.01 1.36 2.36
            k      6.29 6.31 6.83 6.47 5.62 5.22 4.73 4.81 4.56 3.82   <- declining
```

The phase model is essentially sound and merely runs 30-50% small. The offsets model is missing
an additive term — which is consistent with what the code does: the offsets budget
(`make3DOffsets.c:367`) contains only `interpRangeSigma` (the `.sr` band, a *local* neighbourhood
scatter computed after removing a local plane) plus a DEM term and the baseline *parameter*
covariance, and omits `offsets.sigmaRresidual`, the rparams tie-point fit residual. The phase
budget's equivalent term, `min(6*PI, vhParam->sigma)`, is present
(`make3DMosaic.c:639,703`).

### Normalised residual z

```
z = (difference - median(difference)) / sqrt(formal^2 + refSE^2)
```

`refSE = s1Std / sqrt(nValid)` is the reference's own standard error. Including it makes `z`
self-consistent: if *both* error models are right, `z` has unit scale. Normalising by the formal
error alone would charge all of the reference's scatter to the product.

`frac(|z|>3)` should be ~0.3% for calibrated errors. The released offsets give 51-55%.

### Semivariogram

`gamma(h)` of the difference at 1.6-51.2 km lags. This distinguishes a **short-correlation**
error (matching noise, which spatial averaging suppresses) from a **long-correlation** one
(ionosphere or an orbit ramp, which it does not).

The distinction is not academic. The `.sr` band that supplies the offsets' entire random-error
budget is a local neighbourhood statistic computed *after a local plane is removed*
(`Cullst/cullStats.c:78-90`), so it cannot represent long-wavelength error by construction —
and the ionospheric correction applied at `make3DOffsets.c:326,337` carries no uncertainty at
all. If the variogram is still climbing at tens of km, no per-pixel local statistic will ever
capture the missing error.

### Cross-product comparison

`--commonMask` intersects the products' valid footprints and recomputes everything on the
identical pixel set, warning if that removes more than 10% of any product. Without it, a
coverage difference masquerades as an error difference. The pairwise `NISAR_a - NISAR_b` scales
are reported directly.

## The resolution caveat

The release COGs are built with `RESAMPLING=AVERAGE` (`setupquarters.py:1066-1068`), so
`--overviewLevel 2` reads an 8x8-**averaged** NISAR field. A natively-coarse 1600 m mosaic
point-samples instead — `mosaic3d` does not average sub-pixels. These are different quantities:

- **Formal errors are near-invariant** to overview level (averaging 64 spatially smooth `ex`
  values returns roughly the same `ex`).
- **The difference scale is not** — averaging suppresses white noise up to 8x and long-wavelength
  error not at all.

That asymmetry is useful rather than annoying: the ratio between the native-1600 m and
averaged-1600 m difference scales is an independent measure of the same short-vs-long
correlation question the variogram answers.

Use `--overviewLevel 2` to reproduce the published numbers; `-1` (native) for coarse test
mosaics.

## Mask

Built from the reference alone, so the same pixels are used for every product and run:

```
isfinite(s1.vv) & s1.vv <= maxSpeed & s1Std <= maxS1Std & nValid >= minEpochs, eroded `erode` px
```

Defaults (`maxS1Std=inf`, `minEpochs=1`, `erode=0`) reproduce the notebook's `vv <= 50`. Tighten
with `--maxS1Std`/`--minEpochs` to restrict to the most stable ground — but report both, because
changing the mask moves the numbers and comparing a tightened result against the published table
would be misleading.

## Outputs

- `<tag>.stats.yaml` — every statistic, plus provenance: each product's `mosaic-000.000.log`
  flag block (`timeThresh`, `timeThreshPhase`, `pairOverCount`, `Feather length`, ...), the
  `mosaic3d` binary's mtime/size/SHA-256, and the S1 cache identity.
- `<tag>.summary.csv` — one row per (product, component).
- `--writeMaps` — `diff` and `z` GeoTIFFs on the product grid, for finding bad areas in QGIS.

The provenance block matters. `mosaic3d` lives outside a git repo, and these products carry no
record of which error model produced them. An attempt to date a code change from output
statistics alone reached the wrong conclusion twice; see "Dating trap" in
`GIT64/mosaicSource/CLAUDE.md`. Read the log, never infer.
