#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compareToS1.py

Compare NISAR velocity mosaics against the cached Sentinel-1 reference and
report how far mosaic3d's formal errors are from the observed differences.

Build the cache first with makeS1Reference.  See Documents/compareToS1.md for
what each statistic means and how to read the output.

The question this answers is not "how big is the difference" but "is the formal
error a correct description of it".  Three things make that answerable that a
single ratio cannot:

  kGlobal   scale(difference) / median(formal).  1 means calibrated.
  kDecile   the same ratio computed within deciles of the formal error.  Flat
            means one global scale factor fixes it; rising means the formal
            error is not even a good RELATIVE predictor and a per-pixel term is
            missing.
  variogram whether the unexplained difference is short-correlation (matching
            noise, which spatial averaging suppresses) or long-correlation
            (ionosphere/ramp, which it does not).  A local neighbourhood
            statistic such as the .sr band can only ever capture the first.

Products are compared on a COMMON mask built from the reference alone, so the
same pixels are used for every product; --commonMask additionally intersects
the products' own valid footprints, since a coverage difference would otherwise
masquerade as an error difference.

Usage:
    compareToS1 --product ROOT --label NAME [--product ... --label ...]
                [--s1Cache DIR] [--overviewLevel N] [--commonMask]
                [--maxSpeed 50] [--maxS1Std INF] [--minEpochs 1] [--erode 0]
                [--outDir DIR] [--tag NAME] [--writeMaps] [--noFigures]
"""
import argparse
import csv
import hashlib
import os
import re
import sys
from datetime import datetime, timezone

import numpy as np
import yaml

from nisarerrors import velErrorStats as ves
from nisarerrors.makeS1Reference import readReference, readBand, DEFAULT_OUT_DIR

MOSAIC3D_BIN = '/home/ian/bin/x86_64/mosaic3d'
COMPONENTS = ['vx', 'vy']
# mosaic3d log keys worth recording: which error model produced the product.
LOG_KEYS = ['timeThresh', 'timeThreshPhase', 'pairOverCount', 'useSquint',
            'noMask', 'Feather length', 'NoVh Flag', 'rOffset Flag',
            'No3D Flag', 'ThreeD Off Flag', 'TimeOverlapFlag',
            'verticalCorrection', 'inputFile', 'sepAscDesc']


def productPaths(root, band):
    """Expand a release root containing '*' into the band's tif."""
    if '*' in root:
        return root.replace('*', band) + '.tif'
    return f'{root}_{band}.tif'


def readProduct(root, bands, overviewLevel=-1):
    """Read a NISAR release product.

    overviewLevel >= 0 reads that GDAL overview instead of full resolution,
    which is how the published comparison numbers were produced (both sides at
    overview 2).  Note this AVERAGES 8x8 at level 2, so it is not the same
    quantity as a natively-coarse mosaic -- see Documents/compareToS1.md.
    """
    from osgeo import gdal
    gdal.UseExceptions()
    out, gt, wkt = {}, None, None
    for band in bands:
        path = productPaths(root, band)
        ds = gdal.Open(path)
        b = ds.GetRasterBand(1)
        if overviewLevel is not None and overviewLevel >= 0:
            if overviewLevel >= b.GetOverviewCount():
                raise RuntimeError(f'{path} has {b.GetOverviewCount()} '
                                   f'overviews; level {overviewLevel} '
                                   f'requested')
            src = b.GetOverview(overviewLevel)
            fx = ds.RasterXSize / src.XSize
            fy = ds.RasterYSize / src.YSize
            g = ds.GetGeoTransform()
            thisGt = (g[0], g[1] * fx, g[2], g[3], g[4], g[5] * fy)
        else:
            src = b
            thisGt = ds.GetGeoTransform()
        a = src.ReadAsArray().astype(np.float64)
        nd = b.GetNoDataValue()
        if nd is not None:
            a[a == nd] = np.nan
        # mosaic3d's own sentinel, which survives some conversions
        a[a < -1.0e9] = np.nan
        out[band] = a
        if gt is None:
            gt, wkt = thisGt, ds.GetProjection()
    return out, gt, wkt


def gridsMatch(gtA, shapeA, gtB, shapeB, tol=0.01):
    """Same shape and geotransform to within tol metres.

    The tolerance is deliberately tight.  A native 1600 m mosaic and the S1
    overview-2 grid have the SAME shape and the SAME origin but pixel sizes of
    1600 and 1600.211/1600.467 -- a 0.013% difference that accumulates to 0.125
    px in x and 0.5 px in y at the far edges.  A loose tolerance would call
    those equal and subtract them directly, silently misregistering the
    southern rows by a whole pixel.
    """
    if shapeA != shapeB:
        return False
    return all(abs(a - b) <= tol for a, b in zip(gtA, gtB))


def regridToTarget(src, srcGt, targetGt, targetShape, method='nearest'):
    """Sample the REFERENCE onto the product's grid.

    Always this direction: resampling the product would smooth away the very
    per-pixel noise being measured.
    """
    ny, nx = targetShape
    xT = targetGt[0] + (np.arange(nx) + 0.5) * targetGt[1]
    yT = targetGt[3] + (np.arange(ny) + 0.5) * targetGt[5]
    cols = (xT - srcGt[0]) / srcGt[1] - 0.5
    rows = (yT - srcGt[3]) / srcGt[5] - 0.5
    sny, snx = src.shape
    if method == 'nearest':
        ci = np.clip(np.round(cols).astype(int), 0, snx - 1)
        ri = np.clip(np.round(rows).astype(int), 0, sny - 1)
        out = src[np.ix_(ri, ci)].astype(np.float64).copy()
        bad = ((cols < -0.5) | (cols > snx - 0.5))[None, :] | \
              ((rows < -0.5) | (rows > sny - 0.5))[:, None]
        out[bad] = np.nan
        return out
    # bilinear, NaN-propagating
    c0 = np.floor(cols).astype(int)
    r0 = np.floor(rows).astype(int)
    fc = cols - c0
    fr = rows - r0
    out = np.full(targetShape, np.nan)
    valid = ((c0 >= 0) & (c0 < snx - 1))[None, :] & \
            ((r0 >= 0) & (r0 < sny - 1))[:, None]
    c0c = np.clip(c0, 0, snx - 2)
    r0c = np.clip(r0, 0, sny - 2)
    v00 = src[np.ix_(r0c, c0c)]
    v01 = src[np.ix_(r0c, c0c + 1)]
    v10 = src[np.ix_(r0c + 1, c0c)]
    v11 = src[np.ix_(r0c + 1, c0c + 1)]
    w = ((1 - fr)[:, None] * (1 - fc)[None, :] * v00 +
         (1 - fr)[:, None] * fc[None, :] * v01 +
         fr[:, None] * (1 - fc)[None, :] * v10 +
         fr[:, None] * fc[None, :] * v11)
    out[valid] = w[valid]
    return out


def parseMosaicLog(productRoot):
    """Recover the flags that produced a product from its sector log.

    Provenance must be READ, never inferred from the output statistics -- doing
    the latter is what "Dating trap" in mosaicSource/CLAUDE.md records as having
    produced two wrong conclusions.
    """
    relDir = os.path.dirname(productRoot)
    velDir = os.path.dirname(relDir) if relDir.endswith('release') else relDir
    for cand in ('pieces/mosaic-000.000.log', 'pieces/mosaic-002.002.log'):
        path = os.path.join(velDir, cand)
        if os.path.exists(path):
            break
    else:
        return {'logFound': False}
    flags = {'logFound': True, 'logPath': path}
    with open(path, errors='ignore') as fp:
        for line in fp:
            if not line.startswith(';'):
                continue
            m = re.match(r';\s*([^:]+?)\s*:\s*(\S.*?)\s*$', line)
            if m and any(k.lower() in m.group(1).lower() for k in LOG_KEYS):
                flags[m.group(1).strip()] = m.group(2).strip()
    return flags


def binaryProvenance(path=MOSAIC3D_BIN):
    """mosaic3d is not in a git repo, so identity is mtime+size+sha256."""
    if not os.path.exists(path):
        return {'path': path, 'exists': False}
    h = hashlib.sha256()
    with open(path, 'rb') as fp:
        for block in iter(lambda: fp.read(1 << 20), b''):
            h.update(block)
    st = os.stat(path)
    return {'path': path, 'exists': True, 'size': st.st_size,
            'mtime': datetime.fromtimestamp(st.st_mtime,
                                            timezone.utc).isoformat(),
            'sha256': h.hexdigest()}


def analyse(products, cache, args):
    """Run every product against the reference.  Returns the results dict."""
    s1Gt = tuple(cache['geoTransform'])
    results = {}

    # First pass: read products, put the reference on each product's grid.
    prepared = {}
    for label, root in products:
        bands = COMPONENTS + ['vv', 'ex', 'ey']
        data, gt, _ = readProduct(root, bands, overviewLevel=args.overviewLevel)
        shape = data['vx'].shape
        if gridsMatch(s1Gt, cache['mean']['vx'].shape, gt, shape):
            regrid = 'none'
            ref = {k: cache['mean'][k] for k in ('vx', 'vy', 'vv')}
            refStd = {k: cache['std'][k] for k in ('vx', 'vy')}
            nValid = cache['nValid']
        else:
            regrid = args.regridMethod
            ref = {k: regridToTarget(cache['mean'][k], s1Gt, gt, shape, regrid)
                   for k in ('vx', 'vy', 'vv')}
            refStd = {k: regridToTarget(cache['std'][k], s1Gt, gt, shape,
                                        regrid) for k in ('vx', 'vy')}
            nValid = regridToTarget(cache['nValid'], s1Gt, gt, shape, 'nearest')
        mask = ves.stableMask(ref['vv'], s1Std=refStd['vx'], nValid=nValid,
                              maxSpeed=args.maxSpeed, maxS1Std=args.maxS1Std,
                              minEpochs=args.minEpochs, erode=args.erode)
        prepared[label] = dict(root=root, data=data, gt=gt, shape=shape,
                               ref=ref, refStd=refStd, nValid=nValid,
                               mask=mask, regrid=regrid)

    # Optional intersection of the products' own valid footprints.
    if args.commonMask and len(prepared) > 1:
        shapes = {p['shape'] for p in prepared.values()}
        if len(shapes) != 1:
            print('WARNING: products have different shapes; skipping '
                  '--commonMask')
        else:
            inter = None
            ownValid = {}
            for label, p in prepared.items():
                # This product's OWN valid data within the stable mask -- the
                # right denominator.  Measuring the loss against the stable
                # mask instead would report the same figure for every product
                # and say nothing about whether they cover the same ground.
                v = p['mask'] & np.isfinite(p['data']['vx']) & \
                    np.isfinite(p['data']['ex']) & (p['data']['ex'] > 0)
                ownValid[label] = v
                inter = v if inter is None else (inter & v)
            nStable = int(next(iter(prepared.values()))['mask'].sum())
            print(f'stable mask: {nStable} px; intersection: '
                  f'{int(inter.sum())} px')
            for label, p in prepared.items():
                own = int(ownValid[label].sum())
                p['maskOwn'] = own
                p['mask'] = inter
                lost = 1.0 - inter.sum() / max(own, 1)
                p['commonMaskLoss'] = float(lost)
                p['coverageOfStable'] = float(own / max(nStable, 1))
                print(f'  {label}: covers {100*own/max(nStable,1):.1f}% of the '
                      f'stable mask; --commonMask drops a further '
                      f'{100*lost:.1f}% of its own pixels')
                if lost > 0.10:
                    print(f'  WARNING: {label} loses {100*lost:.1f}% to the '
                          f'intersection -- these products are not covering '
                          f'the same ground, so uncorrected comparisons '
                          f'between them are partly a coverage difference')

    pixelKm = abs(prepared[products[0][0]]['gt'][1]) / 1000.0
    for label, _ in products:
        p = prepared[label]
        entry = {'root': p['root'], 'regrid': p['regrid'],
                 'shape': list(p['shape']), 'pixelSizeM': abs(p['gt'][1]),
                 'maskPixels': int(p['mask'].sum()),
                 'flags': parseMosaicLog(p['root'])}
        if 'commonMaskLoss' in p:
            entry['commonMaskLoss'] = p['commonMaskLoss']
            entry['maskPixelsOwn'] = p['maskOwn']
        for comp, err in (('vx', 'ex'), ('vy', 'ey')):
            diff = p['ref'][comp] - p['data'][comp]
            refSE = p['refStd'][comp] / np.sqrt(np.maximum(p['nValid'], 1))
            entry[comp] = ves.compareOne(
                diff, p['data'][err], p['mask'], refSE=refSE,
                pixelSizeKm=pixelKm, withVariogram=not args.noVariogram)
        results[label] = entry

    # Product-vs-product: the direct measure of whether combining helped.
    if len(products) > 1:
        pair = {}
        labels = [l for l, _ in products]
        for i, a in enumerate(labels):
            for b in labels[i + 1:]:
                if prepared[a]['shape'] != prepared[b]['shape']:
                    continue
                m = prepared[a]['mask'] & prepared[b]['mask']
                d = {}
                for comp in COMPONENTS:
                    delta = prepared[a]['data'][comp] - prepared[b]['data'][comp]
                    dv = delta[m & np.isfinite(delta)]
                    if dv.size:
                        d[comp] = {'mad': float(ves.robustScale(dv)),
                                   'std': float(np.std(dv)),
                                   'median': float(np.median(dv))}
                pair[f'{a} - {b}'] = d
        results['_pairwise'] = pair
    return results, prepared


def summaryTable(results, products):
    """The user's table, in its original shape."""
    rows = []
    for label, _ in products:
        r = results[label]
        row = {'product': label, 'n': r['maskPixels']}
        for comp, err in (('vx', 'ex'), ('vy', 'ey')):
            c = r.get(comp, {})
            if 'diff' not in c:
                continue
            row[f'{comp}_diffStd'] = c['diff']['std']
            row[f'{comp}_diffMad'] = c['diff']['mad']
            row[f'{err}_mean'] = c['formal']['mean']
            row[f'{err}_median'] = c['formal']['median']
            row[f'{comp}_k'] = c['kGlobal']
        rows.append(row)
    return rows


def printTable(rows):
    print(f"\n{'product':22s}{'n':>10s}"
          f"{'dStd x':>9s}{'dStd y':>9s}"
          f"{'ex mean':>9s}{'ey mean':>9s}"
          f"{'k x':>8s}{'k y':>8s}")
    for r in rows:
        print(f"{r['product']:22s}{r['n']:10d}"
              f"{r.get('vx_diffStd', np.nan):9.2f}"
              f"{r.get('vy_diffStd', np.nan):9.2f}"
              f"{r.get('ex_mean', np.nan):9.3f}"
              f"{r.get('ey_mean', np.nan):9.3f}"
              f"{r.get('vx_k', np.nan):8.2f}"
              f"{r.get('vy_k', np.nan):8.2f}")


def main():
    p = argparse.ArgumentParser(
        epilog='Part of the nisarErrors package.', description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--product', action='append', default=[],
                   help='Release product root containing "*" (repeatable)')
    p.add_argument('--label', action='append', default=[],
                   help='Label for each --product, in the same order')
    p.add_argument('--s1Cache', default=DEFAULT_OUT_DIR,
                   help=f'S1 reference cache dir (default: {DEFAULT_OUT_DIR})')
    p.add_argument('--overviewLevel', type=int, default=-1,
                   help='GDAL overview level to read products at; -1 = native '
                        '(default). Use 2 to reproduce the notebook numbers.')
    p.add_argument('--maxSpeed', type=float, default=50.0,
                   help='Stable mask: max S1 speed, m/yr (default 50)')
    p.add_argument('--maxS1Std', type=float, default=np.inf,
                   help='Stable mask: max S1 inter-annual sigma (default inf, '
                        'i.e. the notebook-equivalent mask)')
    p.add_argument('--minEpochs', type=int, default=1,
                   help='Stable mask: min contributing S1 epochs (default 1)')
    p.add_argument('--erode', type=int, default=0,
                   help='Stable mask: erode N pixels (default 0)')
    p.add_argument('--commonMask', action='store_true',
                   help='Intersect the products valid footprints too')
    p.add_argument('--regridMethod', default='nearest',
                   choices=['nearest', 'linear'])
    p.add_argument('--noVariogram', action='store_true')
    p.add_argument('--outDir', default='.', help='Where to write outputs')
    p.add_argument('--tag', default='compareToS1', help='Output name stem')
    p.add_argument('--writeMaps', action='store_true',
                   help='Write difference / formal / z GeoTIFFs')
    args = p.parse_args()

    if not args.product:
        p.error('at least one --product is required')
    if len(args.label) != len(args.product):
        args.label = args.label + [f'product{i}' for i in
                                   range(len(args.label), len(args.product))]
    products = list(zip(args.label, args.product))

    cache = readReference(args.s1Cache)
    results, prepared = analyse(products, cache, args)

    rows = summaryTable(results, products)
    printTable(rows)

    if '_pairwise' in results and results['_pairwise']:
        print('\nproduct-vs-product (robust scale of the difference, m/yr):')
        for k, v in results['_pairwise'].items():
            sx = v.get('vx', {}).get('mad', np.nan)
            sy = v.get('vy', {}).get('mad', np.nan)
            print(f'  {k:34s} vx {sx:7.2f}   vy {sy:7.2f}')

    os.makedirs(args.outDir, exist_ok=True)
    meta = {
        'generated': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'args': {k: (None if v is None else
                     (float(v) if isinstance(v, float) and np.isfinite(v)
                      else ('inf' if isinstance(v, float) else v)))
                 for k, v in vars(args).items()},
        's1Cache': {'root': cache['root'],
                    'built': cache['meta'].get('built'),
                    'nEpochs': cache['meta'].get('nEpochs'),
                    'grid': cache['meta'].get('grid')},
        'mosaic3d': binaryProvenance(),
        'results': results,
    }
    statsPath = os.path.join(args.outDir, f'{args.tag}.stats.yaml')
    with open(statsPath, 'w') as fp:
        yaml.safe_dump(meta, fp, sort_keys=False, default_flow_style=False)
    csvPath = os.path.join(args.outDir, f'{args.tag}.summary.csv')
    if rows:
        keys = sorted({k for r in rows for k in r})
        keys = ['product', 'n'] + [k for k in keys if k not in ('product', 'n')]
        with open(csvPath, 'w', newline='') as fp:
            w = csv.DictWriter(fp, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
    print(f'\nwrote {statsPath}')
    print(f'wrote {csvPath}')

    if args.writeMaps:
        from nisarerrors.makeS1Reference import writeGeoTiff
        from osgeo import gdal
        gdal.UseExceptions()
        for label, _ in products:
            pr = prepared[label]
            wkt = gdal.Open(productPaths(pr['root'], 'vx')).GetProjection()
            for comp, err in (('vx', 'ex'), ('vy', 'ey')):
                d = pr['ref'][comp] - pr['data'][comp]
                z = ves.normalisedResidual(
                    d, pr['data'][err],
                    refSE=pr['refStd'][comp] / np.sqrt(
                        np.maximum(pr['nValid'], 1)))
                for arr, nm in ((d, f'diff_{comp}'), (z, f'z_{comp}')):
                    a = np.where(pr['mask'], arr, np.nan)
                    writeGeoTiff(os.path.join(
                        args.outDir, f'{args.tag}.{label}.{nm}.tif'),
                        a, pr['gt'], wkt, -2.0e9)
        print(f'wrote maps to {args.outDir}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
