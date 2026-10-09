#!/usr/bin/env python3
"""
squintCarrierTests: tests with existing NISAR RUNW data of how squint, coregistration and the
ionosphere enter zero-Doppler InSAR phase.

Theory and results: Documents/squintCoregistrationIonosphere.md.  Sub-commands:

  squintDoppler  recompute squint from orbit + Doppler LUT and compare with geolocationGrid (Fig. 1)
  carrierMaps    predicted vs observed carrier phase on rock for one sub-frame (Fig. 2)
  rock           high-pass / band-pass slope of corrected phase on rubbersheet offset, rock (Fig. 3)
  premise        moving-ice test of the squinted-LOS (-useSquint) premise (Fig. 4)
  crossing       direction of phase/phase vs range/range mosaic3d solutions (Fig. 5)

Sub-frames are given as track-N/ORBIT_FRAME relative to --projectDir; each needs H5/<RUNW>.h5,
offsetSims/offsets.geom{,.mask}.vrt, range.offsets.vrt, azimuth.offsets.vrt and geodat26x16.geojson.
"""
import argparse
import glob
import json
import os
from concurrent.futures import ProcessPoolExecutor

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from osgeo import gdal
from scipy.interpolate import RegularGridInterpolator, CubicHermiteSpline
from scipy.ndimage import gaussian_filter
from scipy.optimize import brentq

gdal.UseExceptions()
C = 299792458.0
DEFAULT_PROJECT_DIR = '/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject'


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def loadRunw(frameDir):
    """RUNW layers, rubbersheet pixelOffsets and scalars for one GrIMP sub-frame directory."""
    f = sorted(glob.glob(f'{frameDir}/H5/NISAR*RUNW*.h5'))[0]
    h = h5py.File(f, 'r')
    s = h['science/LSAR/RUNW/swaths/frequencyA']
    ig, po, m = s['interferogram'], s['pixelOffsets'], h['science/LSAR/RUNW/metadata']
    d = dict(file=f, t=ig['zeroDopplerTime'][()], r=ig['slantRange'][()],
             poT=po['zeroDopplerTime'][()], poR=po['slantRange'][()], dem=ig['digitalElevationModel'][()])
    for k in ['unwrappedPhase', 'ionospherePhaseScreen', 'coherenceMagnitude', 'connectedComponents']:
        d[k] = ig['HH/' + k][()]
    at = po['HH/alongTrackOffset'][()].astype(float)
    at[at == 0] = np.nan                                   # exact zeros are fill
    d['alongTrackOffset'] = at                             # applied rubbersheet shift, m, sec - ref
    d['lam'] = C / s['centerFrequency'][()]
    d['fdc'] = float(np.nanmedian(m['processingInformation/parameters/common/frequencyA/dopplerCentroid'][()]))
    d['vg'] = float(np.nanmedian(m['geolocationGrid/groundTrackVelocity'][1]))
    d['passDir'] = h['science/LSAR/identification/orbitPassDirection'][()].decode()
    h.close()
    return d


def loadVrtBand(path, band=1):
    """Band of a SAR-coordinate VRT plus its pixel-centre (time, range) axes."""
    ds = gdal.Open(path)
    gt = ds.GetGeoTransform()
    b = ds.GetRasterBand(band)
    a = b.ReadAsArray().astype(np.float64)
    nd = b.GetNoDataValue()
    a[a < -1e9] = np.nan
    if nd is not None:
        a[a == nd] = np.nan
    nr, nc = a.shape
    return a, gt[3] + (np.arange(nr) + 0.5) * gt[5], gt[0] + (np.arange(nc) + 0.5) * gt[1]


def interpAxes(z, t, r, tq, rq, method='linear'):
    it = RegularGridInterpolator((t, r), z, method=method, bounds_error=False, fill_value=np.nan)
    T, R = np.meshgrid(tq, rq, indexing='ij')
    return it((T, R))


def ncSmooth(z, w, sigma):
    """Normalised (mask-aware) Gaussian smoothing; returns (smoothed, weight)."""
    num = gaussian_filter(np.where(w, z, 0.0), sigma)
    den = gaussian_filter(w.astype(float), sigma)
    with np.errstate(invalid='ignore', divide='ignore'):
        return num / den, den


def slopeFit(x, y):
    g = np.isfinite(x) & np.isfinite(y)
    if g.sum() < 500:
        return np.nan, np.nan
    return float(np.polyfit(x[g], y[g], 1)[0]), float(np.corrcoef(x[g], y[g])[0, 1])


def demeanPerComponent(z, sel, cc, minPix=2000):
    out = np.full(z.shape, np.nan)
    for cid in np.unique(cc[sel]):
        q = sel & (cc == cid)
        if q.sum() >= minPix:
            out[q] = z[q] - np.nanmedian(z[q])
    return out


# ---------------------------------------------------------------------------
# §6.1 squint from the Doppler centroid
# ---------------------------------------------------------------------------

def squintFromDoppler(h5File, heightIndex=1):
    """Per slant-range column at mid-azimuth: geolocationGrid squint, squint recomputed from orbit at
    the native-Doppler time, and the zero-Doppler residual against alongTrackUnitVector (deg)."""
    from pyproj import Transformer
    h = h5py.File(h5File, 'r')
    m = h['science/LSAR/RUNW/metadata']
    g = m['geolocationGrid']
    lam = C / h['science/LSAR/RUNW/swaths/frequencyA/centerFrequency'][()]
    look = +1 if h['science/LSAR/identification/lookDirection'][()].decode().lower() == 'left' else -1
    o = m['orbit/reference']
    pos = CubicHermiteSpline(o['time'][()], o['position'][()], o['velocity'][()])
    vel = pos.derivative()
    fdcRg = np.nanmedian(m['processingInformation/parameters/common/frequencyA/dopplerCentroid'][()], axis=0)
    j = g['zeroDopplerTime'].shape[0] // 2
    t0 = g['zeroDopplerTime'][j]
    sr = g['slantRange'][()]
    epsg = int(g['epsg'][()])
    toEcef = Transformer.from_crs(epsg, 4978, always_xy=True)
    toLL = Transformer.from_crs(epsg, 4326, always_xy=True)
    X, Y = g['coordinateX'][heightIndex, j, :], g['coordinateY'][heightIndex, j, :]
    H = g['heightAboveEllipsoid'][heightIndex]
    heading = lambda e, n: np.degrees(np.arctan2(e, n))
    wrap = lambda a: (a + 180) % 360 - 180
    rows = []
    for k in range(len(sr)):
        T = np.array(toEcef.transform(X[k], Y[k], H))
        lon, lat = toLL.transform(X[k], Y[k])
        la, lo = np.radians(lat), np.radians(lon)
        enu = np.array([[-np.sin(lo), np.cos(lo), 0],
                        [-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)],
                        [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]])
        f = np.interp(k, np.linspace(0, len(sr) - 1, len(fdcRg)), fdcRg)
        doppler = lambda tt: 2 / lam * np.dot(vel(tt), T - pos(tt)) / np.linalg.norm(T - pos(tt))
        tc = brentq(lambda tt: doppler(tt) - f, t0 - 10, t0 + 10)
        tz = brentq(doppler, t0 - 10, t0 + 10)
        losC = enu @ ((pos(tc) - T) / np.linalg.norm(pos(tc) - T))
        losZ = enu @ ((pos(tz) - T) / np.linalg.norm(pos(tz) - T))
        base = heading(g['alongTrackUnitVectorX'][heightIndex, j, k], g['alongTrackUnitVectorY'][heightIndex, j, k])
        sGrid = wrap(heading(g['losUnitVectorX'][heightIndex, j, k], g['losUnitVectorY'][heightIndex, j, k]) - base) - 90 * look
        sOrbit = wrap(heading(losC[0], losC[1]) - base) - 90 * look
        sZero = wrap(heading(losZ[0], losZ[1]) - base) - 90 * look
        rows.append((sr[k], f, sGrid, sOrbit, sZero))
    h.close()
    return np.array(rows)


# ---------------------------------------------------------------------------
# §6.2/§6.3 carrier term on rock
# ---------------------------------------------------------------------------

def rockTest(frameDir, cohMin=0.5, sigHP=50, sigBP=(50, 150)):
    """Slopes (rad/m) of high-pass and band-pass corrected phase, uncorrected phase and screen on
    the rubbersheet along-track offset, over rock."""
    d = loadRunw(frameDir)
    mask, mt, mr = loadVrtBand(f'{frameDir}/offsetSims/offsets.geom.mask.vrt')
    rock = interpAxes(mask, mt, mr, d['t'], d['r'], 'nearest') == 1
    c = interpAxes(d['alongTrackOffset'], d['poT'], d['poR'], d['t'], d['r'])
    cc = d['connectedComponents']
    unw, scr = d['unwrappedPhase'].astype(float), d['ionospherePhaseScreen'].astype(float)
    sel = (cc > 0) & (d['coherenceMagnitude'] > cohMin) & np.isfinite(unw) & np.isfinite(scr) & rock & np.isfinite(c)
    R, Cc = np.meshgrid(np.arange(unw.shape[0]), np.arange(unw.shape[1]), indexing='ij')
    res = dict(frame=frameDir, passDir=d['passDir'], pred=-2 * np.pi * d['fdc'] / d['vg'])
    for name, src in [('corrected', unw - scr), ('uncorrected', unw), ('screen', scr)]:
        y = demeanPerComponent(src, sel, cc)
        ok = np.isfinite(y)
        if ok.sum() < 5000:
            return None
        A = np.c_[np.ones(ok.sum()), R[ok], Cc[ok], d['dem'][ok]]
        p, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
        y[ok] -= A @ p
        x = np.where(ok, c, np.nan)
        hy, hx = y - ncSmooth(y, ok, sigHP)[0], x - ncSmooth(x, ok, sigHP)[0]
        b1y, den = ncSmooth(y, ok, sigBP[0]); b2y, _ = ncSmooth(y, ok, sigBP[1])
        b1x, _ = ncSmooth(x, ok, sigBP[0]); b2x, _ = ncSmooth(x, ok, sigBP[1])
        m2 = ok & (den > 0.3)
        res[name] = dict(highPass=slopeFit(hx[ok], hy[ok]), bandPass=slopeFit((b1x - b2x)[m2][::25], (b1y - b2y)[m2][::25]))
        res['n'] = int(ok.sum())
        res['sdC'] = float(np.nanstd(hx[ok]))
    return res


# ---------------------------------------------------------------------------
# §6.4 moving-ice premise test
# ---------------------------------------------------------------------------

def premiseTest(frameDir, cohMin=0.6, mMin=2.0):
    """Regress D = unw + rangeOffsetPhase - 2*screen on [plane, m - c, c] over moving ice."""
    d = loadRunw(frameDir)
    p = json.load(open(f'{frameDir}/geodat26x16.geojson'))['properties']
    mask, mt, mr = loadVrtBand(f'{frameDir}/offsetSims/offsets.geom.mask.vrt')
    ice = interpAxes(mask, mt, mr, d['t'], d['r'], 'nearest') == 2
    geomA, gt, gr = loadVrtBand(f'{frameDir}/offsetSims/offsets.geom.vrt', 2)
    geomR, _, _ = loadVrtBand(f'{frameDir}/offsetSims/offsets.geom.vrt', 1)
    azOff, _, _ = loadVrtBand(f'{frameDir}/azimuth.offsets.vrt')
    rgOff, _, _ = loadVrtBand(f'{frameDir}/range.offsets.vrt')
    m = interpAxes((azOff - geomA) * p['SLCAzimuthPixelSize'], gt, gr, d['t'], d['r'])
    rOff = interpAxes((rgOff - geomR) * p['SLCRangePixelSize'], gt, gr, d['t'], d['r'])
    c = interpAxes(d['alongTrackOffset'], d['poT'], d['poR'], d['t'], d['r'])
    cc = d['connectedComponents']
    unw, scr = d['unwrappedPhase'].astype(float), d['ionospherePhaseScreen'].astype(float)
    Dfull = unw - 4 * np.pi / d['lam'] * rOff - 2 * scr
    sel = (cc > 0) & (d['coherenceMagnitude'] > cohMin) & np.isfinite(Dfull) & np.isfinite(c) & np.isfinite(m) & ice & (np.abs(m) > mMin)
    D = demeanPerComponent(Dfull, sel, cc, minPix=3000)
    ok = np.isfinite(D)
    if ok.sum() < 20000:
        return None
    med = np.nanmedian(D[ok]); mad = 1.4826 * np.nanmedian(np.abs(D[ok] - med))
    ok &= np.abs(D - med) < 4 * mad
    R, Cc = np.meshgrid(np.arange(unw.shape[0]), np.arange(unw.shape[1]), indexing='ij')
    A = np.c_[np.ones(ok.sum()), R[ok], Cc[ok], (m - c)[ok], c[ok]]
    pp, *_ = np.linalg.lstsq(A, D[ok], rcond=None)
    se = np.sqrt(np.diag(np.linalg.inv(A.T @ A) * np.var(D[ok] - A @ pp) * 25))   # x25: ~5x5 px correlation
    return dict(frame=frameDir, passDir=d['passDir'], n=int(ok.sum()), pred=-2 * np.pi * d['fdc'] / d['vg'],
                coefC=float(pp[4]), seC=float(se[4]), coefMC=float(pp[3]), seMC=float(se[3]),
                c95=float(np.percentile(np.abs(c[ok]), 95)), c=c[ok][::5], D=(D[ok] - A[:, :3] @ pp[:3])[::5])


# ---------------------------------------------------------------------------
# §6.5 crossing-orbit direction
# ---------------------------------------------------------------------------

def readMosaic(base):
    ds = gdal.Open(base + '.vx.tif')
    gt = ds.GetGeoTransform()
    vx = ds.ReadAsArray().astype(float)
    vy = gdal.Open(base + '.vy.tif').ReadAsArray().astype(float)
    ex = gdal.Open(base + '.ex.tif').ReadAsArray().astype(float)
    ey = gdal.Open(base + '.ey.tif').ReadAsArray().astype(float)
    bad = (vx < -1e9) | (vy < -1e9) | ~np.isfinite(vx) | ~np.isfinite(vy) | (ex <= 0) | (ey <= 0)
    z = vx + 1j * vy
    z[bad] = np.nan
    ny, nx = vx.shape
    X, Y = np.meshgrid(gt[0] + (np.arange(nx) + .5) * gt[1], gt[3] + (np.arange(ny) + .5) * gt[5])
    return z, np.hypot(ex, ey), X / 1e3, Y / 1e3


def crossingDirection(runDir, snrMin=10, blockKm=25, nBoot=300, seed=0):
    """Block-median direction of phaseOnly and phaseSquint relative to rangeOnly (deg, ± bootstrap)."""
    P, eP, X, Y = readMosaic(f'{runDir}/phaseOnly')
    S, _, _, _ = readMosaic(f'{runDir}/phaseSquint')
    R, eR, _, _ = readMosaic(f'{runDir}/rangeOnly')
    g = np.isfinite(P) & np.isfinite(S) & np.isfinite(R) & (np.abs(R) / eR > snrMin) & (np.abs(P) / eP > snrMin)
    blk = np.floor(X[g] / blockKm) * 1e5 + np.floor(Y[g] / blockKm)
    ub = np.unique(blk)
    out = dict(n=int(g.sum()), nBlocks=len(ub))
    rng = np.random.default_rng(seed)
    for name, A in [('phase', P), ('squint', S)]:
        dt = np.degrees(np.angle(A[g] * np.conj(R[g])))
        bmed = np.array([np.median(dt[blk == b]) for b in ub])
        nb = np.array([(blk == b).sum() for b in ub])
        w = np.minimum(nb, 200).astype(int)
        boots = [np.median(np.repeat(bmed[k], w[k])) for k in (rng.choice(len(ub), len(ub)) for _ in range(nBoot))]
        out[name] = (float(np.median(np.repeat(bmed, w))), float(np.std(boots)))
    return out


def crossingRunScript(projectDir, ascTrack, descTrack, outDir, dem):
    """Write inputFile and a csh script running the three mosaic3d variants for one crossing."""
    os.makedirs(outDir, exist_ok=True)
    lines = []
    for t in (ascTrack, descTrack):
        for d in sorted(glob.glob(f'{projectDir}/{t}/[0-9]*_0000')):
            need = [f'{d}/phase.uw.interp.vrt', f'{d}/motion/baseline.26x16.yaml',
                    f'{d}/motion/rBaseline.deltabp.yaml', f'{d}/motion/az.est.yaml']
            if all(os.path.exists(n) for n in need):
                lines.append(f'{d}/phase.uw.interp.vrt {d}/geodat26x16.geojson {d}/motion/baseline.26x16.yaml 12 1. '
                             f'{d}/azimuth.offsets {d}/motion/az.est.yaml {d}/range.offsets {d}/motion/rBaseline.deltabp.yaml')
    with open(f'{outDir}/inputFile', 'w') as fp:
        fp.write(f'0 0 0 0 .2 .2\n{len(lines)}\n;\n' + '\n'.join(lines) + '\n')
    common = '-GTiff -ompThreads 6 -noVh -fl 0'
    with open(f'{outDir}/runCrossing.csh', 'w') as fp:
        fp.write('setenv OPENBLAS_NUM_THREADS 1\n'
                 f'mosaic3d {common} -noRangeRows -noAzimuthRows inputFile {dem} phaseOnly >& phaseOnly.out &\n'
                 f'mosaic3d {common} -noRangeRows -noAzimuthRows -useSquint inputFile {dem} phaseSquint >& phaseSquint.out &\n'
                 f'mosaic3d {common} -3dOff -no3d -offsets -noPhaseRows -noAzimuthRows -timeThresh 25 inputFile {dem} rangeOnly >& rangeOnly.out &\n'
                 'wait\n')
    return len(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _mapFrames(func, frames, projectDir, nProc):
    with ProcessPoolExecutor(nProc) as ex:
        return [r for r in ex.map(func, [os.path.join(projectDir, f) for f in frames]) if r is not None]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('command', choices=['squintDoppler', 'rock', 'premise', 'crossingSetup', 'crossing'])
    ap.add_argument('items', nargs='*', help='RUNW h5 files (squintDoppler), sub-frames (rock, premise), '
                    'asc desc tracks (crossingSetup) or mosaic3d run dirs (crossing)')
    ap.add_argument('--projectDir', default=DEFAULT_PROJECT_DIR)
    ap.add_argument('--outDir', default='.')
    ap.add_argument('--dem', default='/Volumes/insar7/ian/gimp/gimp2/270m/dem.gimp2.270m')
    ap.add_argument('--nProc', type=int, default=6)
    a = ap.parse_args()
    if a.command == 'squintDoppler':
        for f in a.items:
            r = squintFromDoppler(f)
            print(f'{os.path.basename(f)[17:40]} max|grid-orbit| {np.max(np.abs(r[:, 2] - r[:, 3])):.4f} deg, '
                  f'zero-Doppler residual {np.median(r[:, 4]):+.3f}, squint vs zero-Doppler LOS {np.median(r[:, 2] - r[:, 4]):.3f}')
    elif a.command == 'rock':
        for r in _mapFrames(rockTest, a.items, a.projectDir, a.nProc):
            print(f"{r['frame'][-18:]} pred {r['pred']:+.2f}  HP corrected {r['corrected']['highPass'][0]:+.2f}  "
                  f"BP corrected {r['corrected']['bandPass'][0]:+.2f}  BP screen {r['screen']['bandPass'][0]:+.2f}")
    elif a.command == 'premise':
        for r in _mapFrames(premiseTest, a.items, a.projectDir, a.nProc):
            print(f"{r['frame'][-18:]} coef c {r['coefC']:+.3f}±{r['seC']:.3f}  coef m-c {r['coefMC']:+.3f}±{r['seMC']:.3f}  "
                  f"(premise +{-r['pred']:.2f})")
    elif a.command == 'crossingSetup':
        n = crossingRunScript(a.projectDir, a.items[0], a.items[1], a.outDir, a.dem)
        print(f'{n} frames; run: cd {a.outDir} && csh runCrossing.csh')
    elif a.command == 'crossing':
        for d in a.items:
            r = crossingDirection(d)
            print(f"{d}: n={r['n']} blocks={r['nBlocks']}  phase {r['phase'][0]:+.2f}±{r['phase'][1]:.2f}  "
                  f"useSquint {r['squint'][0]:+.2f}±{r['squint'][1]:.2f} deg")


if __name__ == '__main__':
    main()
