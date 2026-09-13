"""Validate a velocity mosaic's formal errors against the Sentinel-1 reference.

The acceptance test used to calibrate the crossing-pair over-count parameter
(see mosaicSource/Documents/mosaic3d.md).  Generalised so any mosaic can be
checked the same way -- in particular a new solver such as the per-pixel
normal-equation crossing-orbit method in
mosaicSource/Documents/crossingOrbitRedundancy.md.

Definitions, which matter:

  e        per-pixel formal error, the ex/ey bands
  d        S1 - NISAR difference at a pixel, vx and vy pooled
  RMS(e)   sqrt(mean(e^2)) -- the POPULATION prediction from the formal errors.
           Variances average, standard deviations do not, so mean(e) and
           median(e) are both wrong here and bias k high by 25-40%.
  std(d)   observed standard deviation of the differences
  k        std(d) / RMS(e).  k = 1 is calibrated; k > 1 means the reported
           errors are too small.

All statistics are taken over stable ground (S1 speed <= maxSpeed) on the
footprint common to every product being compared, so numbers from one run are
directly comparable with another.

Usage:

    python -m nisarerrors.validateMosaic TAG [TAG ...]
    python -m nisarerrors.validateMosaic --coverage TAG

where TAG is either a full path to a merged product prefix or a bare name that
is resolved against --mergedDir.  A product prefix P is expected to have
P_vx.tif, P_vy.tif, P_ex.tif, P_ey.tif alongside it.
"""
import argparse
import numpy as np

from nisarerrors.makeS1Reference import readReference, DEFAULT_OUT_DIR
from nisarerrors.plotErrCal import loadPair

DEFAULT_MERGED = ('/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject/'
                  'Release/velocity/errCal/merged')

LEVELS = [0.683, 0.900, 0.950, 0.990, 0.999]
GAUSS = {0.683: 1.000, 0.900: 1.645, 0.950: 1.960, 0.990: 2.576, 0.999: 3.291}


def rmsE(e):
    """Population prediction from the per-pixel formal errors."""
    return float(np.sqrt(np.nanmean(e ** 2)))


def loadProducts(tags, mergedDir, refDir):
    """Load each tag paired with the S1 reference. Returns {tag: (nisar, ref)}."""
    cache = readReference(refDir)
    out = {}
    for tag in tags:
        path = tag if '/' in tag else f'{mergedDir}/{tag}'
        d, r = loadPair(path, cache)[:2]
        out[tag] = (d, r)
    return out


def commonMask(loaded, maxSpeed):
    """Stable-ground pixels valid in EVERY product, so runs are comparable."""
    common = None
    for d, r in loaded.values():
        v = (np.isfinite(r['vv']) & (r['vv'] <= maxSpeed) &
             np.isfinite(d['vx']) & np.isfinite(d['ex']) & (d['ex'] > 0))
        common = v if common is None else (common & v)
    return common


def calibration(loaded, mask):
    """RMS(e), std(d) and k for each product on the shared mask."""
    rows = {}
    for tag, (d, r) in loaded.items():
        e = np.concatenate([d['ex'][mask], d['ey'][mask]])
        dif = np.concatenate([(r['vx'] - d['vx'])[mask],
                              (r['vy'] - d['vy'])[mask]])
        f, sd = rmsE(e), float(np.nanstd(dif))
        rows[tag] = (f, sd, sd / f)
    return rows


def coverage(loaded, mask):
    """Coverage and confidence multipliers, with the overall scale divided out.

    z is rescaled so RMS(z) = 1 exactly -- i.e. assuming the scale has already
    been calibrated -- which isolates the SHAPE of the error distribution from
    the scale question.  Departures from the Gaussian row are then purely about
    the tail, not about whether the sigmas are the right size.
    """
    rows = {}
    for tag, (d, r) in loaded.items():
        e = np.concatenate([d['ex'][mask], d['ey'][mask]])
        dif = np.concatenate([(r['vx'] - d['vx'])[mask],
                              (r['vy'] - d['vy'])[mask]])
        z = dif / e
        z = z / np.sqrt(np.nanmean(z ** 2))
        a = np.abs(z[np.isfinite(z)])
        rows[tag] = ([float((a < m).mean()) for m in (1, 2, 3)],
                     [float(np.quantile(a, L)) for L in LEVELS])
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('tags', nargs='+', help='merged product prefixes or full paths')
    p.add_argument('--mergedDir', default=DEFAULT_MERGED)
    p.add_argument('--refDir', default=DEFAULT_OUT_DIR)
    p.add_argument('--maxSpeed', type=float, default=50.0,
                   help='stable ground cutoff, S1 m/yr (default 50)')
    p.add_argument('--coverage', action='store_true',
                   help='also report coverage and confidence multipliers')
    args = p.parse_args()

    loaded = loadProducts(args.tags, args.mergedDir, args.refDir)
    mask = commonMask(loaded, args.maxSpeed)
    print(f'common footprint {int(mask.sum())} px, S1 speed <= {args.maxSpeed} m/yr\n')

    print(f'{"product":>28s} {"RMS(e)":>9s} {"std(d)":>9s} {"k":>7s}')
    for tag, (f, sd, k) in calibration(loaded, mask).items():
        print(f'{tag:>28s} {f:9.3f} {sd:9.3f} {k:7.2f}')

    if args.coverage:
        print(f'\n{"product":>28s} {"<1s":>7s} {"<2s":>7s} {"<3s":>7s}   ' +
              ' '.join(f'{L*100:>6.1f}%' for L in LEVELS))
        for tag, (cov, mult) in coverage(loaded, mask).items():
            print(f'{tag:>28s} ' + ' '.join(f'{c*100:6.1f}%' for c in cov) +
                  '   ' + ' '.join(f'{m:6.2f}' for m in mult))
        print(f'{"Gaussian":>28s} ' +
              ' '.join(f'{v:6.1f}%' for v in (68.3, 95.4, 99.7)) +
              '   ' + ' '.join(f'{GAUSS[L]:6.2f}' for L in LEVELS))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
