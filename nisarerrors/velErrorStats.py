#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
velErrorStats.py

Statistics for comparing a NISAR velocity mosaic against a reference velocity
field, used to calibrate mosaic3d's formal errors (see compareToS1 and
Documents/compareToS1.md).

Deliberately numpy/xarray only -- no nisardev, no grimpfunc, no network -- so
it stays importable from a notebook without EarthData credentials.  The driver
does the I/O; this module does the arithmetic.

Two methodological points drive the design.

1) Comparing a STD of the difference against a MEAN of the formal error is not
   like-for-like: the first is outlier-sensitive, the second is not.  Every
   scale here is therefore reported both ways (std and 1.4826*MAD), and every
   scale is reported both raw and after removing the median difference -- a
   one-signed bias is not something a formal error bar is meant to reproduce,
   so folding it into the comparison overstates the error.

2) A single ratio cannot distinguish "the sigmas are globally under-scaled" from
   "the sigmas are fine except in localised bad areas".  kDecile answers that:
   if k is flat across deciles of the formal error, one global factor fixes it;
   if k rises, the formal error is not even a good RELATIVE predictor and a
   per-pixel term is missing.  Those imply different fixes.
"""
import numpy as np

# Median-absolute-deviation to Gaussian-sigma conversion.
MAD_TO_SIGMA = 1.4826
DEFAULT_PERCENTILES = [1, 5, 16, 25, 50, 75, 84, 95, 99]
DEFAULT_LAGS_KM = [1.6, 3.2, 6.4, 12.8, 25.6, 51.2]


def robustScale(x):
    """1.4826 * MAD -- a Gaussian-consistent scale insensitive to outliers."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    return MAD_TO_SIGMA * np.median(np.abs(x - np.median(x)))


def spreadFromPercentiles(x, lo=16, hi=84):
    """(p84 - p16)/2 -- a second robust scale, less sensitive than MAD to
    a sharply peaked core."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    a, b = np.percentile(x, [lo, hi])
    return (b - a) / 2.0


def scaleSummary(x, percentiles=DEFAULT_PERCENTILES):
    """Location and scale of a sample, robust and non-robust side by side."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {'n': 0}
    med = float(np.median(x))
    out = {
        'n': int(x.size),
        'mean': float(np.mean(x)),
        'median': med,
        'std': float(np.std(x)),
        'rmse': float(np.sqrt(np.mean(x ** 2))),
        'mad': float(robustScale(x)),
        'p84p16': float(spreadFromPercentiles(x)),
    }
    # Scales after removing the median: what remains once a common-mode bias,
    # which no formal error bar models, is taken out.
    xc = x - med
    out['stdDemedianed'] = float(np.std(xc))
    out['madDemedianed'] = float(robustScale(xc))
    out['pct'] = {str(p): float(v)
                  for p, v in zip(percentiles, np.percentile(x, percentiles))}
    return out


def stableMask(s1Speed, s1Std=None, nValid=None, maxSpeed=50.0,
               maxS1Std=np.inf, minEpochs=1, erode=0):
    """Pixels stable enough to treat the reference as truth.

    Built from the REFERENCE only, never from the product under test, so the
    same pixels are used for every run and every product and results stay
    comparable.

    erode removes `erode` pixels from the mask boundary, suppressing edge
    pixels where the two fields' footprints disagree.
    """
    m = np.isfinite(s1Speed) & (s1Speed <= maxSpeed)
    if s1Std is not None and np.isfinite(maxS1Std):
        m &= np.isfinite(s1Std) & (s1Std <= maxS1Std)
    if nValid is not None and minEpochs > 1:
        m &= np.isfinite(nValid) & (nValid >= minEpochs)
    for _ in range(int(erode)):
        m = _erodeOnce(m)
    return m


def _erodeOnce(m):
    """4-connected binary erosion, without pulling in scipy."""
    e = m.copy()
    e[1:, :] &= m[:-1, :]
    e[:-1, :] &= m[1:, :]
    e[:, 1:] &= m[:, :-1]
    e[:, :-1] &= m[:, 1:]
    return e


def normalisedResidual(diff, formal, refSE=None, deMedian=True):
    """z = (d - median(d)) / sqrt(formal^2 + refSE^2).

    Including the reference's own standard error in the denominator makes z
    self-consistent: if BOTH error models are right, z has unit scale.
    Normalising by the formal error alone would charge all of the reference's
    scatter to the product under test.
    """
    d = np.asarray(diff, dtype=float)
    if deMedian:
        good = np.isfinite(d)
        if good.any():
            d = d - np.median(d[good])
    denom = np.asarray(formal, dtype=float) ** 2
    if refSE is not None:
        denom = denom + np.asarray(refSE, dtype=float) ** 2
    denom = np.sqrt(denom)
    with np.errstate(divide='ignore', invalid='ignore'):
        z = np.where(denom > 0, d / denom, np.nan)
    return z


def kFactor(diff, formal, useRobust=True):
    """How wrong the sigmas are: scale(d) / median(formal).

    k = 1 means calibrated; k = 11 means the formal errors are 11x too small.
    """
    d = np.asarray(diff, dtype=float)
    f = np.asarray(formal, dtype=float)
    good = np.isfinite(d) & np.isfinite(f) & (f > 0)
    if good.sum() < 10:
        return np.nan
    scale = robustScale(d[good]) if useRobust else np.std(d[good])
    return float(scale / np.median(f[good]))


def kByDecile(diff, formal, nBins=10, useRobust=True):
    """k computed within deciles of the formal error.

    THE decision-relevant statistic.  Flat across deciles => the formal error is
    a good relative predictor and one global scale factor fixes it.  Rising =>
    it is not, and a per-pixel term is missing rather than a constant.
    """
    d = np.asarray(diff, dtype=float)
    f = np.asarray(formal, dtype=float)
    good = np.isfinite(d) & np.isfinite(f) & (f > 0)
    d, f = d[good], f[good]
    if d.size < nBins * 10:
        return []
    edges = np.percentile(f, np.linspace(0, 100, nBins + 1))
    out = []
    for i in range(nBins):
        lo, hi = edges[i], edges[i + 1]
        sel = (f >= lo) & (f <= hi) if i == nBins - 1 else (f >= lo) & (f < hi)
        if sel.sum() < 10:
            out.append({'bin': i, 'n': int(sel.sum())})
            continue
        scale = robustScale(d[sel]) if useRobust else float(np.std(d[sel]))
        out.append({'bin': i, 'n': int(sel.sum()),
                    'formalLo': float(lo), 'formalHi': float(hi),
                    'formalMedian': float(np.median(f[sel])),
                    'diffScale': float(scale),
                    'k': float(scale / np.median(f[sel]))})
    return out


def semivariogram(field, mask, pixelSizeKm, lagsKm=None, maxPairs=200000,
                  rng=None):
    """Empirical semivariogram of `field` over `mask`.

    gamma(h) = 0.5 * mean( (v(x) - v(x+h))^2 ) for separations near h.

    This is what distinguishes a short-correlation error (white matching noise,
    which spatial averaging suppresses) from a long-correlation one (an
    ionospheric or baseline ramp, which averaging does not touch).  A variogram
    that saturates by a few km says the missing error is local; one still
    climbing at tens of km says it is not, and no per-pixel local statistic
    such as a neighbourhood scatter can ever capture it.

    Sampled by random pixel pairs at each lag rather than exhaustively.
    """
    if lagsKm is None:
        lagsKm = DEFAULT_LAGS_KM
    if rng is None:
        rng = np.random.default_rng(12345)   # fixed: results reproducible
    v = np.asarray(field, dtype=float)
    ok = mask & np.isfinite(v)
    rows, cols = np.nonzero(ok)
    if rows.size < 100:
        return []
    out = []
    for lagKm in lagsKm:
        lagPix = max(1, int(round(lagKm / pixelSizeKm)))
        n = min(maxPairs, rows.size)
        idx = rng.choice(rows.size, size=n, replace=False)
        r0, c0 = rows[idx], cols[idx]
        # random axis and sign, so the estimate is direction-averaged
        ang = rng.uniform(0, 2 * np.pi, size=n)
        r1 = r0 + np.round(lagPix * np.sin(ang)).astype(int)
        c1 = c0 + np.round(lagPix * np.cos(ang)).astype(int)
        inside = ((r1 >= 0) & (r1 < v.shape[0]) &
                  (c1 >= 0) & (c1 < v.shape[1]))
        r0, c0, r1, c1 = r0[inside], c0[inside], r1[inside], c1[inside]
        if r0.size < 50:
            out.append({'lagKm': float(lagKm), 'nPairs': 0})
            continue
        both = ok[r1, c1]
        r0, c0, r1, c1 = r0[both], c0[both], r1[both], c1[both]
        if r0.size < 50:
            out.append({'lagKm': float(lagKm), 'nPairs': 0})
            continue
        dv = v[r0, c0] - v[r1, c1]
        # The classical estimator is a MEAN of squared differences, so a
        # handful of extreme pixels dominate it -- on point-sampled coarse
        # mosaics this makes gamma(h) non-monotonic and unusable.  Report the
        # robust scale alongside it and prefer that one for interpretation.
        out.append({'lagKm': float(lagKm), 'lagPix': int(lagPix),
                    'nPairs': int(dv.size),
                    'gamma': float(0.5 * np.mean(dv ** 2)),
                    'sigmaAtLag': float(np.sqrt(np.mean(dv ** 2))),
                    'sigmaAtLagRobust': float(robustScale(dv))})
    return out


def compareOne(diff, formal, mask, refSE=None, pixelSizeKm=None,
               lagsKm=None, withVariogram=True):
    """Full statistic set for one (product, component).

    diff   : reference - product
    formal : the product's formal error for that component
    refSE  : the reference's own standard error (sigma/sqrt(N)), or None
    """
    d = np.asarray(diff, dtype=float)
    f = np.asarray(formal, dtype=float)
    use = mask & np.isfinite(d) & np.isfinite(f) & (f > 0)
    if use.sum() < 10:
        return {'n': int(use.sum())}

    se = None if refSE is None else np.asarray(refSE, dtype=float)
    out = {
        'n': int(use.sum()),
        'diff': scaleSummary(d[use]),
        'formal': scaleSummary(f[use]),
        'kGlobal': kFactor(d[use], f[use]),
        'kGlobalStd': kFactor(d[use], f[use], useRobust=False),
        'kDecile': kByDecile(d[use], f[use]),
    }
    if se is not None:
        seUse = se[use]
        seUse = seUse[np.isfinite(seUse)]
        if seUse.size:
            out['refSE'] = {'median': float(np.median(seUse)),
                            'mean': float(np.mean(seUse))}
            # Reference contribution removed in quadrature.  Reported with the
            # fraction it represents, because when it is a large part of the
            # total this is a small difference of large numbers and the value
            # should not be quoted as if it were well determined.
            dScale = robustScale(d[use])
            seMed = float(np.median(seUse))
            resid = dScale ** 2 - seMed ** 2
            out['refSE']['fractionOfDiffVar'] = float(
                min(1.0, seMed ** 2 / dScale ** 2)) if dScale > 0 else np.nan
            out['productScaleRefRemoved'] = (float(np.sqrt(resid))
                                             if resid > 0 else 0.0)

    z = normalisedResidual(d, f, refSE=se)
    zUse = z[use]
    zUse = zUse[np.isfinite(zUse)]
    if zUse.size:
        out['z'] = {
            'mad': float(robustScale(zUse)),
            'std': float(np.std(zUse)),
            'median': float(np.median(zUse)),
            'fracGt1': float(np.mean(np.abs(zUse) > 1)),
            'fracGt2': float(np.mean(np.abs(zUse) > 2)),
            'fracGt3': float(np.mean(np.abs(zUse) > 3)),
            'pct': {str(p): float(v) for p, v in
                    zip(DEFAULT_PERCENTILES,
                        np.percentile(zUse, DEFAULT_PERCENTILES))},
        }
    if withVariogram and pixelSizeKm is not None:
        out['variogram'] = semivariogram(d, use, pixelSizeKm, lagsKm=lagsKm)
    return out


def selfTest(n=400, sigma=2.0, seed=7):
    """Inject d ~ N(0, sigma) with formal == sigma and confirm the statistics
    recover a calibrated result.  If this fails the machinery is wrong, and any
    conclusion drawn from real data with it would be too."""
    rng = np.random.default_rng(seed)
    d = rng.normal(0.0, sigma, size=(n, n))
    f = np.full_like(d, sigma)
    mask = np.ones_like(d, dtype=bool)
    r = compareOne(d, f, mask, pixelSizeKm=1.6, withVariogram=False)
    checks = {
        'kGlobal~1': abs(r['kGlobal'] - 1.0) < 0.05,
        'zMad~1': abs(r['z']['mad'] - 1.0) < 0.05,
        'fracGt2~0.045': abs(r['z']['fracGt2'] - 0.0455) < 0.01,
        'diffMad~sigma': abs(r['diff']['mad'] - sigma) < 0.05 * sigma,
    }
    for k, v in checks.items():
        print(f'  {k}: {"ok" if v else "FAIL"}')
    return all(checks.values())


if __name__ == '__main__':
    import sys
    print('velErrorStats self-test:')
    sys.exit(0 if selfTest() else 1)
