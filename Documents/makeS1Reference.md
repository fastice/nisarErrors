# makeS1Reference

Builds and caches a Sentinel-1 reference velocity field for calibrating `mosaic3d`'s formal
errors. Used by [compareToS1](compareToS1.md).

```bash
makeS1Reference                 # build (skips if the cache exists)
makeS1Reference --force         # rebuild
makeS1Reference --verify        # rebuild and check against the cache
```

The reference is the temporal mean of the NSIDC-0725 annual Greenland velocity mosaics. Building
it means re-reading the whole remote stack over the network, so it is cached once.

## What is cached, and why three things

| files | source | why |
|---|---|---|
| `*_mean_{vx,vy,vv}.tif` | `vFull.mean()` | the reference velocity |
| `*_std_{vx,vy,vv}.tif` | `vFull.stdev()` | S1's usable uncertainty indicator + the stability criterion |
| `*_nValid.tif` | `vFull.numberValid()` | contributing epochs per pixel |
| `*.yaml` | — | provenance sidecar |

**The stdev, not S1's own formal errors.** The NSIDC-0725 products predate the current
processing updates, so their error bands are not a meaningful reference; the notebook this is
lifted from reads them with `useErrors=False` for the same reason. The year-to-year scatter is
the usable indicator, and independent comparison against ~50 GPS points in northern Greenland
puts S1 at ~1 m/yr or better over stable ground. That is small against the differences being
measured (up to ~10 m/yr), so the reference does not need a modelled error budget — only to be
visible so it is never mistaken for the product's error.

**`nValid` is not optional.** `mean(skipna=True)` will return a one-epoch mean where only one
epoch is valid; those pixels have to be excluded from the stable mask. `sqrt(nValid)` is also
what turns the stdev into a standard error of the mean.

## Grid

Cached on the grid it was *read* on — `overviewLevel=2`, 948 x 1712 — with no resampling at build
time. That is the canonical S1 representation, and it makes a comparison against a 200 m NISAR
product read at the same overview level exactly resampling-free. Putting the reference onto any
other grid is `compareToS1`'s job.

Note the overview grid is **not** a clean 1600 m: its implied pixel size is
`200 * 7585/948 = 1600.211` by `200 * 13700/1712 = 1600.467`, because neither raster dimension
divides by 8. It shares the release grid's NW corner at `(-659100, -639100)` and diverges from a
clean 1600 m grid by up to 0.125 px in x and 0.5 px in y at the far edges. `compareToS1` handles
that with an explicit nearest-neighbour regrid of the reference onto the product grid.

## Provenance sidecar

Records the granule URLs, epoch list, grid, package versions, build time and per-file SHA-256.
`--verify` rebuilds and asserts agreement, and also re-checks the hashes to catch a cache
modified after the fact.

This is not bureaucracy. `mosaic3d` is not in a git repo and its products carried no record of
which error model produced them; an attempt to date a code change from output statistics alone
reached the wrong conclusion twice before a direct rerun settled it. See "Dating trap" in
`GIT64/mosaicSource/CLAUDE.md`.

## Implementation notes worth knowing

**GDAL, not rioxarray, for the writes.** `rioxarray 0.13.1`'s `to_raster` passes a `DataArray` to
xarray's `encode_cf_variable`, which now requires a `Variable`, so every `rio.to_raster` call
raises `AttributeError: 'DataArray' object has no attribute '_data'` in this environment. That
disables `nisardev.writeCloudOptGeo` generally, not just here. Writing through GDAL sidesteps it
and keeps the georeferencing explicit; upgrading rioxarray underneath `nisardev` would have a far
wider blast radius.

**`vv` is written, never derived on read-back.** A mean-of-speeds is not the speed-of-means, so
recomputing `sqrt(vx^2+vy^2)` would silently substitute a different quantity.

**Pixel centres vs corners.** `mosaic.template.yaml` and xarray coordinates use pixel *centres*;
GeoTIFF georeferencing is corner-based. `geoTransformOf()` subtracts the half pixel. Getting this
wrong shifts everything by 800 m at 1600 m spacing, which would look exactly like an
error-magnitude problem at comparison time.

**The login check is advisory.** `NASALogin.check_cookie()` fetches
`urs.earthdata.nasa.gov/profile` and has been observed to return `False` against a cookie jar
that serves NSIDC COGs perfectly well. Gating on it would block a working setup, so it only
warns; a genuinely bad credential fails the first read with a better message. The interactive
login path needs a GUI and is not attempted — refresh the cookie from the notebook if a read
does fail.

**`numberValid()` returns a wrapped `nisarVelSeries`**, not the bare `DataArray` its docstring
describes — it carries the same `_applyInTime` decorator as `mean`/`stdev`. `bandArray()` accepts
either form.
