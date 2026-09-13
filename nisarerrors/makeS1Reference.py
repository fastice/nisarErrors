#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
makeS1Reference.py

Build and cache a Sentinel-1 reference velocity map for calibrating mosaic3d's
formal errors (see compareToS1 in this package, and
Documents/makeS1Reference.md).

The reference is the temporal MEAN of the NSIDC-0725 annual Greenland velocity
mosaics.  Recomputing it means re-reading the whole remote stack, so it is built
once and cached to local COGs.

Three things are cached, not one:

  mean   vx, vy, vv  -- the reference velocity itself.
  stdev  vx, vy, vv  -- the year-to-year scatter.  This is the usable indicator
                        of S1's uncertainty: the NSIDC-0725 products predate the
                        current processing updates, so THEIR OWN formal error
                        bands are not a meaningful reference (the notebook this
                        is lifted from reads them with useErrors=False for the
                        same reason).  The inter-annual sigma also defines which
                        ground is stable enough to compare over.
  nValid             -- per-pixel count of contributing annual epochs.  Needed
                        because mean(skipna=True) will happily return a one-epoch
                        mean; those pixels must be excluded from the stable mask,
                        and sqrt(nValid) is what turns the stdev into a standard
                        error of the mean.

The cache is written on the grid it was READ on (overviewLevel=2 by default,
948 x 1712), with no resampling at build time -- that is the canonical S1
representation, and it keeps a comparison against a 200 m NISAR product read at
the same overview level exactly resampling-free.  Putting S1 onto some other
grid is compareToS1's job.

Every build writes a YAML sidecar recording the granule URLs, epochs, grid,
package versions and per-file SHA-256.  Without it the cache is an undated blob
whose provenance has to be reverse-engineered from its own statistics -- a
failure mode documented under "Dating trap" in mosaicSource/CLAUDE.md.

Usage:
    makeS1Reference [--outDir DIR] [--overviewLevel N] [--firstDate D]
                    [--lastDate D] [--product P] [--verify] [--force]
"""
import argparse
import getpass
import hashlib
import os
import sys
from datetime import datetime, timezone

import numpy as np
import yaml

DEFAULT_OUT_DIR = ('/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject/'
                   'Notebooks/calValReport/data/Greenland/S1reference')
DEFAULT_PRODUCT = 'NSIDC-0725'
DEFAULT_FIRST_DATE = '2014-01-01'
DEFAULT_LAST_DATE = '2026-01-01'
DEFAULT_OVERVIEW_LEVEL = 2
COOKIE_FILE = '.grimp_download_cookiejar.txt'
VEL_BANDS = ['vx', 'vy', 'vv']
# Counts are small positive integers; -1 is unambiguous and survives float32.
NVALID_NODATA = -1.0


def cacheRoot(product, date1, date2, overviewLevel):
    """Root name for the cached files, e.g.
    S1ref_NSIDC-0725_2014-12-01_2023-11-30_ovr2."""
    return f'S1ref_{product}_{date1}_{date2}_ovr{overviewLevel}'


def sidecarPath(outDir, root):
    return os.path.join(outDir, f'{root}.yaml')


def setupGdalEnv(cookieFile=COOKIE_FILE):
    """Point GDAL's HTTP driver at the EarthData cookie jar.

    Must happen before any remote read.  Mirrors the notebook's setup cell.
    """
    cookiePath = os.path.expanduser(f'~/{cookieFile}')
    os.environ['GDAL_HTTP_COOKIEFILE'] = cookiePath
    os.environ['GDAL_HTTP_COOKIEJAR'] = cookiePath
    return cookiePath


def checkLogin(quiet=False):
    """Report on the EarthData cookie, but do NOT gate on it.

    check_cookie() fetches urs.earthdata.nasa.gov/profile, which can fail for
    reasons unrelated to whether the cookie will authorise a data read -- an
    existing jar has been observed to serve NSIDC COGs perfectly well while
    this check returns False.  Gating on it would block a working setup, so it
    is advisory only: if credentials really are bad, the first remote read
    fails loudly and with a more informative error than we could produce here.

    NASALogin's interactive login path needs a GUI/prompt and cannot run
    headless, so it is not attempted; refresh the cookie from the notebook if
    a read does fail.
    """
    import grimpfunc as grimp
    # grimpfunc/__init__.py re-exports the CLASS, not the module.
    try:
        ok = grimp.NASALogin().check_cookie()
    except Exception as e:
        ok = False
        if not quiet:
            print(f'EarthData cookie check errored ({type(e).__name__}); '
                  f'continuing anyway')
    if not quiet:
        print(f'EarthData cookie check: '
              f'{"valid" if ok else "inconclusive -- continuing"}')
    return ok


def loadS1Series(product=DEFAULT_PRODUCT, firstDate=DEFAULT_FIRST_DATE,
                 lastDate=DEFAULT_LAST_DATE,
                 overviewLevel=DEFAULT_OVERVIEW_LEVEL, quiet=False):
    """Search CMR and load the annual velocity stack.

    Lifted from S1-NISARcomparison.ipynb so the cache is provably the same
    object the published comparison numbers were computed from.

    Returns (vFull, urls).
    """
    import nisardev
    import grimpfunc as grimp

    myUrls = grimp.cmrUrls(mode='nisar')
    myUrls.initialSearch(product=product, firstDate=firstDate,
                         lastDate=lastDate)
    # mode='nisar' forces productGroups to ['vv'], so getCogs returns one URL
    # template per epoch with 'vv' replaced by '*'; readSeriesFromTiff then
    # substitutes each band name in turn.
    cogs = myUrls.getCogs(replace='vv', removeTiff=True)
    if len(cogs) == 0:
        raise RuntimeError(f'CMR search returned no {product} granules for '
                           f'{firstDate}..{lastDate}')
    if not quiet:
        print(f'found {len(cogs)} {product} epochs')

    vFull = nisardev.nisarVelSeries()
    # readSpeed=False -> only vx/vy are fetched and vv is derived per epoch,
    # which is what the published comparison used.
    vFull.readSeriesFromTiff(cogs, url=True, readSpeed=False, useErrors=False,
                             useStack=True, overviewLevel=overviewLevel)
    vFull.loadRemote()
    if not quiet:
        print(f'loaded stack shape {vFull.subset.shape} '
              f'bands {list(vFull.subset.band.data)}')
    return vFull, cogs


def buildReference(vFull):
    """Reduce the stack to (vMean, vSigma, nValid).

    vMean/vSigma are nisarVelSeries instances (nisarBase2D's _applyInTime
    decorator rewraps them); nValid is a bare DataArray of dims (band, y, x).
    """
    vMean = vFull.mean()
    vSigma = vFull.stdev()
    nValid = vFull.numberValid()
    return vMean, vSigma, nValid


def _squeezeTime(obj):
    """Drop the length-1 time dim left by mean()/stdev()/numberValid().

    These are wrapped by nisarBase2D's _applyInTime decorator, so they return
    nisarVelSeries instances rather than bare DataArrays -- but numberValid's
    docstring and body disagree about that, so accept either form.
    """
    xr = getattr(obj, 'subset', obj)
    if 'time' in getattr(xr, 'dims', ()):
        xr = xr.isel(time=0)
    return xr


def bandArray(obj, band):
    """Float array for one band of a velocity object or DataArray."""
    return np.asarray(_squeezeTime(obj).sel(band=band).data, dtype=float)


def sha256Of(path, blockSize=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as fp:
        for block in iter(lambda: fp.read(blockSize), b''):
            h.update(block)
    return h.hexdigest()


def geoTransformOf(xr):
    """GDAL geotransform from an xarray whose x/y coords are pixel CENTRES.

    GeoTIFF georeferencing is corner-based, so half a pixel is subtracted.
    Getting this wrong shifts the whole grid by 800 m at 1600 m spacing, which
    would look exactly like an error-magnitude problem at comparison time.
    """
    x, y = xr.x.data, xr.y.data
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])          # negative when y runs north -> south
    return (float(x[0]) - dx / 2., dx, 0.0,
            float(y[0]) - dy / 2., 0.0, dy)


def writeGeoTiff(path, array, geoTransform, wkt, noData):
    """Write a single-band float32 GeoTIFF via GDAL.

    Deliberately not rioxarray: rioxarray 0.13.1's to_raster hands a DataArray
    to xarray's encode_cf_variable, which now requires a Variable, so every
    rio.to_raster call raises AttributeError in this environment.  Writing
    through GDAL sidesteps that and keeps the georeferencing explicit.
    """
    from osgeo import gdal
    gdal.UseExceptions()
    a = np.asarray(array, dtype=np.float32)
    a = np.where(np.isfinite(a), a, noData).astype(np.float32)
    drv = gdal.GetDriverByName('GTiff')
    ds = drv.Create(path, a.shape[1], a.shape[0], 1, gdal.GDT_Float32,
                    options=['COMPRESS=DEFLATE', 'TILED=YES',
                             'BIGTIFF=IF_SAFER'])
    ds.SetGeoTransform(geoTransform)
    ds.SetProjection(wkt)
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(float(noData))
    band.WriteArray(a)
    band.FlushCache()
    ds.BuildOverviews('AVERAGE', [2, 4, 8])
    ds = None
    return path


def writeReference(vMean, vSigma, nValid, outDir, root, quiet=False):
    """Write mean/stdev/nValid GeoTIFFs.  Returns the list of files written."""
    from osgeo import osr
    os.makedirs(outDir, exist_ok=True)
    written = []

    sq = _squeezeTime(vMean)
    gt = geoTransformOf(sq)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(str(sq.rio.crs).split(':')[-1]))
    wkt = srs.ExportToWkt()

    for velSeries, tag in ((vMean, 'mean'), (vSigma, 'std')):
        xrSq = _squeezeTime(velSeries)
        # vv MUST be written, not derived on read-back: a mean-of-speeds is not
        # the speed-of-means, and recomputing sqrt(vx^2+vy^2) would silently
        # substitute a different quantity.
        for band in VEL_BANDS:
            path = os.path.join(outDir, f'{root}_{tag}_{band}.tif')
            writeGeoTiff(path, xrSq.sel(band=band).data, gt, wkt,
                         getattr(velSeries, 'noDataDict', {}).get(band,
                                                                  -2.0e9))
            written.append(path)

    # Counts: vx and vy are identical in practice; keep vx's.
    nValidPath = os.path.join(outDir, f'{root}_nValid.tif')
    writeGeoTiff(nValidPath, bandArray(nValid, 'vx'), gt, wkt, NVALID_NODATA)
    written.append(nValidPath)

    if not quiet:
        for f in written:
            print(f'  wrote {os.path.basename(f)}')
    return written


def gridDescription(velSeries):
    """Grid of the cache, recorded so a mismatch is caught rather than
    silently interpolated away at comparison time."""
    xr = velSeries.subset
    x, y = xr.x.data, xr.y.data
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    return {'shape': [int(y.size), int(x.size)],
            'dx': dx, 'dy': dy,
            # pixel CENTRES, matching mosaic.template.yaml's convention
            'xMin': float(x.min()), 'xMax': float(x.max()),
            'yMin': float(y.min()), 'yMax': float(y.max()),
            # pixel CORNER, matching the GeoTIFF geotransform origin
            'originCorner': [float(x.min() - dx / 2.),
                             float(y.max() - dy / 2.)],
            'crs': str(xr.rio.crs)}


def writeSidecar(outDir, root, product, urls, vFull, vMean, overviewLevel,
                 written, args, quiet=False):
    import nisardev
    import grimpfunc

    times = vFull.subset.time.data
    epochs = [str(np.datetime_as_string(t, unit='D')) for t in times]

    meta = {
        'product': product,
        'nEpochs': int(len(epochs)),
        'epochs': epochs,
        'date1': str(np.datetime_as_string(vFull.subset.time1.data[0],
                                           unit='D')),
        'date2': str(np.datetime_as_string(vFull.subset.time2.data[-1],
                                           unit='D')),
        'overviewLevel': int(overviewLevel),
        'searchFirstDate': args.firstDate,
        'searchLastDate': args.lastDate,
        'grid': gridDescription(vMean),
        'urls': list(urls),
        'files': {os.path.basename(f): sha256Of(f) for f in written},
        'built': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'builtBy': getpass.getuser(),
        'versions': {
            'nisardev': getattr(nisardev, '__version__', 'unknown'),
            'grimpfunc': getattr(grimpfunc, '__version__', 'unknown'),
            'python': sys.version.split()[0],
        },
    }
    path = sidecarPath(outDir, root)
    with open(path, 'w') as fp:
        yaml.safe_dump(meta, fp, sort_keys=False, default_flow_style=False)
    if not quiet:
        print(f'  wrote {os.path.basename(path)}')
    return path, meta


def readBand(path):
    """Read one GeoTIFF band as a float array with nodata -> NaN.

    Returns (array, geoTransform, wkt).
    """
    from osgeo import gdal
    gdal.UseExceptions()
    ds = gdal.Open(path)
    b = ds.GetRasterBand(1)
    a = b.ReadAsArray().astype(np.float64)
    nd = b.GetNoDataValue()
    if nd is not None:
        a[a == nd] = np.nan
    a[~np.isfinite(a)] = np.nan
    return a, ds.GetGeoTransform(), ds.GetProjection()


def readReference(outDir, root=None, quiet=False):
    """Read a cached reference back.  Used by compareToS1.

    Read with GDAL rather than nisardev so the caller gets plain arrays plus an
    explicit geotransform -- alignment against a product is then something to
    assert, not something to hope for.

    Returns dict with mean/std band arrays, nValid, geoTransform, wkt, meta.
    """
    if root is None:
        cands = [f for f in os.listdir(outDir) if f.endswith('.yaml')]
        if len(cands) != 1:
            raise RuntimeError(f'expected exactly one sidecar in {outDir}, '
                               f'found {len(cands)} -- pass root explicitly')
        root = cands[0].replace('.yaml', '')
    with open(sidecarPath(outDir, root)) as fp:
        meta = yaml.safe_load(fp)

    out = {'root': root, 'meta': meta, 'mean': {}, 'std': {}}
    gt = wkt = None
    for tag in ('mean', 'std'):
        for band in VEL_BANDS:
            a, g, w = readBand(os.path.join(outDir,
                                            f'{root}_{tag}_{band}.tif'))
            out[tag][band] = a
            if gt is None:
                gt, wkt = g, w
            elif g != gt:
                raise RuntimeError(f'{root}_{tag}_{band}.tif geotransform '
                                   f'differs from the rest of the cache')
    out['nValid'], _, _ = readBand(os.path.join(outDir,
                                                f'{root}_nValid.tif'))
    out['geoTransform'], out['wkt'] = gt, wkt
    if not quiet:
        print(f'read cache {root} grid {out["mean"]["vx"].shape} '
              f'dx={gt[1]:.4f}')
    return out


def verifyReference(outDir, root, args, quiet=False):
    """Rebuild over a small window and check it matches the cache.

    Guards against a silently corrupted or stale cache, and against a nisardev
    change altering mean()/stdev() semantics under us.
    """
    cache = readReference(outDir, root, quiet=quiet)
    meta = cache['meta']
    print(f'cache built {meta["built"]} from {meta["nEpochs"]} epochs')

    print('re-reading source stack for verification ...')
    vFull, _ = loadS1Series(product=meta['product'],
                            firstDate=meta['searchFirstDate'],
                            lastDate=meta['searchLastDate'],
                            overviewLevel=meta['overviewLevel'], quiet=quiet)
    rMean, rSigma, rValid = buildReference(vFull)

    ok = True
    for tag, rebuilt in (('mean', rMean), ('std', rSigma)):
        for band in VEL_BANDS:
            a = cache[tag][band]
            b = bandArray(rebuilt, band)
            both = np.isfinite(a) & np.isfinite(b)
            if not both.any():
                print(f'  {tag}.{band}: NO overlapping valid pixels')
                ok = False
                continue
            dMax = float(np.nanmax(np.abs(a[both] - b[both])))
            flag = 'ok' if dMax < 1e-4 else 'FAIL'
            print(f'  {tag}.{band}: max|diff| = {dMax:.3e}  [{flag}]')
            ok &= dMax < 1e-4

    a = cache['nValid']
    b = bandArray(rValid, 'vx')
    both = np.isfinite(a) & np.isfinite(b)
    nBad = int((a[both] != b[both]).sum())
    print(f'  nValid: {nBad} mismatched pixels  '
          f'[{"ok" if nBad == 0 else "FAIL"}]')
    ok &= nBad == 0

    for name, digest in meta['files'].items():
        path = os.path.join(outDir, name)
        actual = sha256Of(path)
        if actual != digest:
            print(f'  {name}: SHA-256 MISMATCH -- cache modified since build')
            ok = False
    return ok


def main():
    p = argparse.ArgumentParser(
        epilog='Part of the nisarErrors package.', description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--outDir', default=DEFAULT_OUT_DIR,
                   help=f'Cache directory (default: {DEFAULT_OUT_DIR})')
    p.add_argument('--product', default=DEFAULT_PRODUCT,
                   help=f'NSIDC product (default: {DEFAULT_PRODUCT})')
    p.add_argument('--firstDate', default=DEFAULT_FIRST_DATE,
                   help=f'Search start (default: {DEFAULT_FIRST_DATE})')
    p.add_argument('--lastDate', default=DEFAULT_LAST_DATE,
                   help=f'Search end (default: {DEFAULT_LAST_DATE})')
    p.add_argument('--overviewLevel', type=int,
                   default=DEFAULT_OVERVIEW_LEVEL,
                   help=f'Overview level to read and cache at '
                        f'(default: {DEFAULT_OVERVIEW_LEVEL} -> 1600 m)')
    p.add_argument('--verify', action='store_true',
                   help='Verify an existing cache against a fresh rebuild')
    p.add_argument('--force', action='store_true',
                   help='Rebuild even if the cache already exists')
    p.add_argument('--quiet', action='store_true')
    args = p.parse_args()

    setupGdalEnv()
    checkLogin(quiet=args.quiet)

    if args.verify:
        cands = [f for f in os.listdir(args.outDir) if f.endswith('.yaml')]
        if not cands:
            print(f'ERROR: no cache found in {args.outDir}')
            return 1
        ok = verifyReference(args.outDir, cands[0].replace('.yaml', ''),
                             args, quiet=args.quiet)
        print('VERIFY PASSED' if ok else 'VERIFY FAILED')
        return 0 if ok else 1

    vFull, urls = loadS1Series(product=args.product, firstDate=args.firstDate,
                               lastDate=args.lastDate,
                               overviewLevel=args.overviewLevel,
                               quiet=args.quiet)
    vMean, vSigma, nValid = buildReference(vFull)

    date1 = str(np.datetime_as_string(vFull.subset.time1.data[0], unit='D'))
    date2 = str(np.datetime_as_string(vFull.subset.time2.data[-1], unit='D'))
    root = cacheRoot(args.product, date1, date2, args.overviewLevel)

    if os.path.exists(sidecarPath(args.outDir, root)) and not args.force:
        print(f'cache {root} already exists in {args.outDir}; '
              f'use --force to rebuild')
        return 0

    written = writeReference(vMean, vSigma, nValid, args.outDir, root,
                             quiet=args.quiet)
    _, meta = writeSidecar(args.outDir, root, args.product, urls, vFull,
                           vMean, args.overviewLevel, written, args,
                           quiet=args.quiet)
    g = meta['grid']
    print(f'\ncached {root}')
    print(f'  {meta["nEpochs"]} epochs {meta["date1"]}..{meta["date2"]}')
    print(f'  grid {g["shape"][1]} x {g["shape"][0]}  '
          f'dx={g["dx"]:.4f} dy={g["dy"]:.4f}  corner={g["originCorner"]}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
