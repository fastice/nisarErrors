#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plotErrCal.py

Reproduce the two calValReport comparison figures for an arbitrary NISAR
velocity mosaic, so error-budget variants can be compared like for like:

  Figure 3 style -- a 3x3 map panel:
      row 1  |S1 - NISAR| for vx, vy, vv, with the stable-region sigma in the title
      row 2  the NISAR formal errors ex, ey, and the propagated speed error ev
      row 3  the L2-SCI-667 slow-flow requirement 1 + 3%|v|
  Figure 4 style -- a GNSS bar chart at the NIU interior sites:
      per site, S1 mean minus GNSS, NISAR minus GNSS, and the NISAR formal error,
      with the +/-(1 + 3%|obs|) tolerance marked.

Follows the notebook S1-NISARcomparison.ipynb (cells 29 and 32) so the output is
directly comparable to the published figures, but reads the cached S1 reference
(makeS1Reference) and plain GeoTIFFs instead of the remote stack, so it needs no
credentials and can be pointed at any number of test mosaics.

Note ev is the SPEED error propagated from ex/ey, not an independent band:
    ev = sqrt(vx^2 ex^2 + vy^2 ey^2) / vv        (vv guarded away from 0)

Usage:
    plotErrCal --product ROOT --label NAME [--product ... --label ...]
               [--s1Cache DIR] [--gpsYaml FILE] [--outDir DIR]
               [--overviewLevel N] [--maxSpeed 50]
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from nisarerrors.makeS1Reference import readReference, DEFAULT_OUT_DIR
from nisarerrors.compareToS1 import readProduct, gridsMatch, regridToTarget
from nisarerrors import velErrorStats as ves

DEFAULT_GPS_YAML = ('/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject/'
                    'Notebooks/calValReport/data/GPS.interior.2526.yaml')
VMAX = 10.0          # shared colour scale, as in the published figure
CMAP = 'Reds'


def speedError(vx, vy, ex, ey, vv, floor=1.0):
    """Propagated speed error.  vv is guarded so slow pixels do not blow up."""
    denom = np.where(np.abs(vv) > floor, vv, np.nan)
    return np.sqrt(vx ** 2 * ex ** 2 + vy ** 2 * ey ** 2) / denom


def loadPair(root, cache, overviewLevel=-1, regridMethod='nearest'):
    """Read a product and put the S1 reference on its grid."""
    bands = ['vx', 'vy', 'vv', 'ex', 'ey']
    data, gt, _ = readProduct(root, bands, overviewLevel=overviewLevel)
    shape = data['vx'].shape
    s1Gt = tuple(cache['geoTransform'])
    if gridsMatch(s1Gt, cache['mean']['vx'].shape, gt, shape):
        ref = {k: cache['mean'][k] for k in ('vx', 'vy', 'vv')}
        refStd = {k: cache['std'][k] for k in ('vx', 'vy')}
        nValid = cache['nValid']
    else:
        ref = {k: regridToTarget(cache['mean'][k], s1Gt, gt, shape,
                                 regridMethod) for k in ('vx', 'vy', 'vv')}
        refStd = {k: regridToTarget(cache['std'][k], s1Gt, gt, shape,
                                    regridMethod) for k in ('vx', 'vy')}
        nValid = regridToTarget(cache['nValid'], s1Gt, gt, shape, 'nearest')
    return data, ref, refStd, nValid, gt


def figureMaps(data, ref, label, outPath, maxSpeed=50.0, mask=None):
    """Figure 3 style 3x3 panel.

    `mask` lets the caller impose a footprint common to every product being
    compared.  Without it each product is scored on its own coverage, and a
    product that additionally fills gaps where the others have nothing is
    penalised for the extra ground rather than credited for it -- the gap-fill
    pixels are the hardest in the scene and dominate any whole-map sigma.
    """
    slow = (np.isfinite(ref['vv']) & (ref['vv'] <= maxSpeed)) if mask is None \
        else mask
    vvSafe = np.where(np.abs(ref['vv']) > 1, ref['vv'], np.nan)
    evF = speedError(data['vx'], data['vy'], data['ex'], data['ey'],
                     data['vv'])
    reqVx = 1 + 0.03 * np.abs(ref['vx'])
    reqVy = 1 + 0.03 * np.abs(ref['vy'])
    reqVv = np.sqrt(ref['vx'] ** 2 * reqVx ** 2 +
                    ref['vy'] ** 2 * reqVy ** 2) / vvSafe

    fig, axes = plt.subplots(3, 3, figsize=(21, 21))
    plt.subplots_adjust(hspace=0.06, wspace=0.02, left=0.01, right=0.99,
                        top=0.95, bottom=0.05)
    fig.suptitle(f'{label}', fontsize=22, y=0.975)

    # Every panel is masked to the stable region, not just the statistics.
    # Outside it the S1-NISAR difference is dominated by REAL velocity change
    # along the coast -- a decade separates the two epochs -- which has nothing
    # to do with the error budget and would otherwise dominate the figure while
    # being excluded from every number quoted on it.
    def show(ax, arr, title):
        im = ax.imshow(np.where(slow, arr, np.nan), vmin=0, vmax=VMAX,
                       cmap=CMAP, origin='upper')
        ax.set_title(title, fontsize=15)
        ax.axis('off')
        return im

    for ax, band in zip(axes[0], ['vx', 'vy', 'vv']):
        d = ref[band] - data[band]
        ok = slow & np.isfinite(d)
        # robust and classical side by side: these distributions have heavy
        # tails (std/MAD ~ 2.5-3), so quoting only sigma overstates the typical
        # disagreement
        sig = np.nanstd(d[ok])
        mad = ves.robustScale(d[ok])
        im = show(ax, np.abs(d),
                  f'|{band}: S1 - NISAR|   sigma={sig:.1f}  MAD={mad:.1f} m/yr')
    for ax, f, nm in zip(axes[1], [data['ex'], data['ey'], evF],
                         ['ex', 'ey', 'ev']):
        ok = slow & np.isfinite(f)
        im = show(ax, f, f'NISAR {nm}   mean={np.nanmean(f[ok]):.1f}  '
                         f'median={np.nanmedian(f[ok]):.1f} m/yr')
    for ax, r, nm in zip(axes[2], [reqVx, reqVy, reqVv], ['vx', 'vy', 'vv']):
        ok = slow & np.isfinite(r)
        im = show(ax, r, f'requirement 1+3%|{nm}|   '
                         f'mean={np.nanmean(r[ok]):.1f} m/yr')
    cax = fig.add_axes([0.25, 0.025, 0.5, 0.012])
    fig.colorbar(im, cax=cax, orientation='horizontal',
                 extend='max', label='m/yr')
    fig.savefig(outPath, dpi=90, bbox_inches='tight')
    plt.close(fig)


def loadGps(path):
    import yaml
    d = yaml.safe_load(open(path))
    names = list(d)
    g = {k: np.array([d[s]['total'][k] for s in names])
         for k in ('x', 'y', 'vx', 'vy', 'speed')}
    g['names'] = names
    return g


def sampleAt(arr, gt, xs, ys):
    """Nearest-pixel sample, NaN outside."""
    out = np.full(len(xs), np.nan)
    ny, nx = arr.shape
    for i, (X, Y) in enumerate(zip(xs, ys)):
        c = int(round((X - gt[0]) / gt[1] - 0.5))
        r = int(round((Y - gt[3]) / gt[5] - 0.5))
        if 0 <= r < ny and 0 <= c < nx:
            out[i] = arr[r, c]
    return out


def figureGps(data, ref, gt, gps, label, outPath):
    """Figure 4 style GNSS bar chart."""
    xs, ys = gps['x'], gps['y']
    evF = speedError(data['vx'], data['vy'], data['ex'], data['ey'],
                     data['vv'])
    n = {b: sampleAt(data[b], gt, xs, ys) for b in ('vx', 'vy', 'vv',
                                                    'ex', 'ey')}
    n['ev'] = sampleAt(evF, gt, xs, ys)
    s = {b: sampleAt(ref[b], gt, xs, ys) for b in ('vx', 'vy', 'vv')}

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    x = np.arange(len(gps['names']))
    w = 0.25
    rows = [('vx', s['vx'] - gps['vx'], n['vx'] - gps['vx'], n['ex'],
             1 + 0.03 * np.abs(gps['vx'])),
            ('vy', s['vy'] - gps['vy'], n['vy'] - gps['vy'], n['ey'],
             1 + 0.03 * np.abs(gps['vy'])),
            ('vv', s['vv'] - gps['speed'], n['vv'] - gps['speed'], n['ev'],
             1 + 0.03 * np.abs(gps['speed']))]
    for ax, (comp, s1d, nd, fe, tol) in zip(axes, rows):
        ax.bar(x - w, s1d, w, color='steelblue', edgecolor='k', alpha=.8,
               label='S1 mean - GPS')
        ax.bar(x, nd, w, color='darkorange', edgecolor='k', alpha=.8,
               label='NISAR - GPS')
        ax.bar(x + w, fe, w, color='seagreen', edgecolor='k', alpha=.8,
               label='NISAR formal error')
        for xi, t in zip(x, tol):
            ax.plot([xi - 1.5 * w, xi + 1.5 * w], [t, t], color='red', lw=1.4)
            ax.plot([xi - 1.5 * w, xi + 1.5 * w], [-t, -t], color='red', lw=1.4)
        ax.axhline(0, color='k', lw=.8, ls='--')
        mae = np.nanmean(np.abs(nd))
        maeS1 = np.nanmean(np.abs(s1d))
        ax.set_title(f'{comp}: MAE S1={maeS1:.2f}  NISAR={mae:.2f} m/yr',
                     fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels(gps['names'], rotation=45)
        ax.set_ylabel(f'delta {comp} (product - GPS) m/yr')
    axes[2].legend(fontsize=8)
    fig.suptitle(label, fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(outPath, dpi=110, bbox_inches='tight')
    plt.close(fig)
    return {'names': gps['names'],
            'nisar': {c: r[2] for c, r in zip(('vx', 'vy', 'vv'),
                                              [(0, 0, rows[0][2]),
                                               (0, 0, rows[1][2]),
                                               (0, 0, rows[2][2])])},
            'mae': {c: float(np.nanmean(np.abs(r[2])))
                    for c, r in zip(('vx', 'vy', 'vv'), rows)},
            'formal': {c: float(np.nanmean(r[3]))
                       for c, r in zip(('vx', 'vy', 'vv'), rows)}}


def main():
    p = argparse.ArgumentParser(
        epilog='Part of the nisarErrors package.', description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--product', action='append', default=[])
    p.add_argument('--label', action='append', default=[])
    p.add_argument('--s1Cache', default=DEFAULT_OUT_DIR)
    p.add_argument('--gpsYaml', default=DEFAULT_GPS_YAML)
    p.add_argument('--outDir', default='figures')
    p.add_argument('--overviewLevel', type=int, default=-1)
    p.add_argument('--maxSpeed', type=float, default=50.0)
    p.add_argument('--commonMask', action='store_true',
                   help='Restrict every figure to the footprint where ALL '
                        'products have data, so they are compared like for like')
    args = p.parse_args()
    if not args.product:
        p.error('--product required')
    if len(args.label) != len(args.product):
        args.label += [f'p{i}' for i in range(len(args.label),
                                              len(args.product))]
    os.makedirs(args.outDir, exist_ok=True)
    cache = readReference(args.s1Cache)
    gps = loadGps(args.gpsYaml)

    # Pass 1: load everything, so a footprint common to all can be built.
    loaded = {}
    for label, root in zip(args.label, args.product):
        loaded[label] = loadPair(root, cache, args.overviewLevel)

    common = None
    if args.commonMask:
        for label, (data, ref, _, _, _) in loaded.items():
            v = (np.isfinite(ref['vv']) & (ref['vv'] <= args.maxSpeed) &
                 np.isfinite(data['vx']) & np.isfinite(data['ex']) &
                 (data['ex'] > 0))
            common = v if common is None else (common & v)
        stable = np.isfinite(list(loaded.values())[0][1]['vv']) & \
            (list(loaded.values())[0][1]['vv'] <= args.maxSpeed)
        print(f'common footprint: {int(common.sum())} px of '
              f'{int(stable.sum())} stable ({100*common.sum()/stable.sum():.1f}%)')
        for label, (data, ref, _, _, _) in loaded.items():
            own = int((np.isfinite(data['vx']) &
                       np.isfinite(ref['vv']) &
                       (ref['vv'] <= args.maxSpeed)).sum())
            print(f'   {label}: own coverage {own} px, '
                  f'{100*(1 - common.sum()/max(own,1)):.1f}% dropped')

    # Pass 2: plot everything on the same pixels.
    summary = {}
    for label, root in zip(args.label, args.product):
        data, ref, refStd, nValid, gt = loaded[label]
        m = os.path.join(args.outDir, f'{label}.maps.png')
        g = os.path.join(args.outDir, f'{label}.gps.png')
        figureMaps(data, ref, label, m, maxSpeed=args.maxSpeed, mask=common)
        summary[label] = figureGps(data, ref, gt, gps, label, g)
        print(f'{label}: wrote {os.path.basename(m)}, {os.path.basename(g)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
