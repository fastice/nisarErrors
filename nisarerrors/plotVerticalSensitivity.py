#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plotVerticalSensitivity.py

Illustrates, for the crossing-orbit InSAR-phase velocity inversion in mosaic3d
(make3DMosaic.c -- computeA/computeVxy, documented in
mosaicSource/Documents/mosaic3d.md "Error Analysis: Crossing-Geometry
Sensitivity", and in Documents/plotVerticalSensitivity.md in this package),
how two distinct error sources propagate into the horizontal velocity
components (vx, vy) as a function of pixel position, for three real
ascending/descending crossing pairs at northern, central, and southern
Greenland latitudes:

  1. Uncompensated vertical motion (vz, e.g. unmodeled tide or SMB
     submergence/emergence) -- the SAME physical v_z^SMB contributes to both
     LOS measurements, scaled by each image's own incidence angle:
     delta_pA = v_z^SMB*cot(psi_A), delta_pD = v_z^SMB*cot(psi_D) (NOT simply
     v_z^SMB*1 -- that only holds at psi=45deg). Plotted as dvx, dvy per
     1 m/yr of uncompensated vz, at each region's real incidence angles.

  2. Independent range/phase measurement error (dr) in EACH LOS measurement
     separately -- NOT common-mode. Plotted as sigma_vx, sigma_vy (m/yr) per
     1 metre of raw, independent, equal range-offset noise in each of the two
     LOS measurements, for a 12-day repeat cycle at each region's real
     incidence angle: sigma_vx/vy = rowNorm(A) * 365.25/(12*sin(psi)), where
     rowNorm(A) is the unit-sigma_p sensitivity (row-norm of the A matrix) and
     the scale factor converts 1 metre of raw 12-day range-offset noise into
     sigma_p (m/yr).

Both use the flat-terrain approximation (B = 0, i.e. D = A) since the
question here is the pure crossing-geometry sensitivity, not DEM-slope
effects.

Geometry for the three example crossing pairs: real ascending/descending
track headings, computed from real orbital state vectors in each image's
geodat26x16.geojson (<projectDir>/track-N/orbit_NNNN/), via
headingAndIncidence() -- a faithful port of common/computeHeading.c's
lookDir-aware cross-track-heading formula (NOT simply the flight-direction/
footprint-long-axis angle, which is a materially different quantity -- see
headingAndIncidence() docstring). PassType ('ascending'/'descending') and
LookDirection ('left' here) are read directly from each geodat, not assumed.

Default example geometry (--projectDir overrides the root; the per-region
track/frame relative paths are GIT64/NISAR-Greenland-specific and not
currently overridable from the command line -- edit REGIONS to point at a
different track/frame for a different project):

    Region   Asc geodat              Desc geodat
    North    track-58/3271_37        track-97/3656_53
    Central  track-131/3517_0001     track-11/3397_0000
    South    track-45/3431_0000      track-126/3512_0002

(H_A, H_D, alpha, incidence angle, and the real gamma marker are all computed
at runtime from these files -- printed at the start of main() for inspection.)

Usage:
    plotVerticalSensitivity [--projectDir DIR] [-o OUTPUT.png] [--show]
"""
import argparse
import os
import numpy as np
from osgeo import osr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from utilities import geodatrxa

DEFAULT_PROJECT_DIR = '/Volumes/insar1/ian/NISAR/realNISAR/newGreenlandProject'

# Real crossing pairs: (region label, asc geodat rel path, desc geodat rel path),
# relative to --projectDir.
REGIONS = [
    ('Northern Greenland', 'track-58/3271_37/geodat26x16.geojson',
     'track-97/3656_53/geodat26x16.geojson'),
    ('Central Greenland', 'track-131/3517_0001/geodat26x16.geojson',
     'track-11/3397_0000/geodat26x16.geojson'),
    ('Southern Greenland', 'track-45/3431_0000/geodat26x16.geojson',
     'track-126/3512_0002/geodat26x16.geojson'),
]

REPEAT_DAYS = 12.0  # NISAR-like repeat cycle, for the dr-error y-axis scaling


def computeA(alpha, beta):
    """2x2 geometric conversion matrix, matching common/initRoutines.c computeA()."""
    invSin2 = 1.0 / np.sin(alpha) ** 2
    cosAlphaBeta = np.cos(alpha + beta)
    cosAlpha = np.cos(alpha)
    cosBeta = np.cos(beta)
    sinBeta = np.sin(beta)
    sinAlphaBeta = np.sin(alpha + beta)
    A = np.empty((2, 2))
    A[0, 0] = invSin2 * (cosBeta - cosAlpha * cosAlphaBeta)
    A[0, 1] = invSin2 * (cosAlphaBeta - cosAlpha * cosBeta)
    A[1, 0] = invSin2 * (sinBeta - sinAlphaBeta * cosAlpha)
    A[1, 1] = invSin2 * (sinAlphaBeta - sinBeta * cosAlpha)
    return A


def headingAndIncidence(geojsonPath):
    """Cross-track heading (deg from north, CW) at scene center, faithfully
    porting common/computeHeading.c's lookDir-aware formula onto real orbital
    state vectors -- via geodatrxa.llzPtToRA(), a verified line-for-line port
    of llToImageNew.c's zero-Doppler Newton solve (llToImageNew.c itself has
    no lookDir conditionals; lookDir only enters in the final atan2 step
    below, exactly as in computeHeading.c).

    NOTE: this is the CROSS-TRACK heading computeA() actually needs, not the
    along-track flight-direction heading -- these differ by a curved-Earth,
    incidence-angle-dependent amount (confirmed ~6-16 deg for this dataset),
    not a clean +/-90, so flight heading cannot be substituted here.

    Returns (heading_deg, incidence_deg, lat, lon).
    """
    g = geodatrxa()
    g.readGeojson(geojsonPath)
    lat, lon = g.corners[4, :]  # CenterLatLon
    dlat = -0.0025 if g.lookdir == 'left' else 0.0025
    r1, a1, t1 = g.llzPtToRA(lat - dlat, lon, 0.0)
    r2, a2, t2 = g.llzPtToRA(lat + dlat, lon, 0.0, initT=t1)
    da = (a2 - a1) * g.slpAz
    Re = g.earthRadLatm(np.radians(lat))
    ReH1, ReH2 = g.ReH(t1), g.ReH(t2)
    R1 = r1 * g.slpRg + g.rNearSLP
    R2 = r2 * g.slpRg + g.rNearSLP
    rho1 = np.arccos((R1 ** 2 - ReH1 ** 2 - Re ** 2) / (-2.0 * Re * ReH1))
    rho2 = np.arccos((R2 ** 2 - ReH2 ** 2 - Re ** 2) / (-2.0 * Re * ReH2))
    dgr = Re * rho2 - Re * rho1
    hAngle = np.arctan2(da, dgr) if g.lookdir == 'right' else np.arctan2(da, -dgr)
    return np.degrees(hAngle) % 360, g.phic, lat, lon


def latLonToPS3413(lat, lon):
    src = osr.SpatialReference(); src.ImportFromEPSG(4326)
    tgt = osr.SpatialReference(); tgt.ImportFromEPSG(3413)
    src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    xform = osr.CoordinateTransformation(src, tgt)
    x, y, _ = xform.TransformPoint(lon, lat)
    return x, y


def sensitivities(alpha_deg, gamma_deg, psi_A_deg, psi_D_deg):
    """For heading difference alpha and pixel-azimuth-relative-to-mean-heading
    gamma (both degrees), return (dvx, dvy, rowNormX, rowNormY):
      dvx, dvy   -- common-mode (vz) response per 1 m/yr uncompensated vz.
                    A common-mode v_z^SMB contributes delta_pA =
                    v_z^SMB*cot(psi_A) and delta_pD = v_z^SMB*cot(psi_D)
                    (not simply v_z^SMB*1, which would only hold at
                    psi=45deg) -- see Documents/plotVerticalSensitivity.md
                    for the derivation.
      rowNormX,Y -- independent-error (dr) response per unit sigma_pA=sigma_pD
                    (unscaled; caller applies the m/yr-per-metre-of-dr factor)
    beta = gamma - alpha/2 (since gamma = alpha/2 + beta by definition).
    """
    alpha = np.radians(alpha_deg)
    beta = np.radians(gamma_deg) - alpha / 2.0
    A = computeA(alpha, beta)
    cotPsiA = 1.0 / np.tan(np.radians(psi_A_deg))
    cotPsiD = 1.0 / np.tan(np.radians(psi_D_deg))
    dvx, dvy = A @ np.array([cotPsiA, cotPsiD])
    rowNormX = np.hypot(A[0, 0], A[0, 1])
    rowNormY = np.hypot(A[1, 0], A[1, 1])
    return dvx, dvy, rowNormX, rowNormY


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--projectDir', default=DEFAULT_PROJECT_DIR,
                  help=f'Root directory containing the track-N/ example geodats '
                       f'(default: {DEFAULT_PROJECT_DIR})')
    p.add_argument('-o', '--output', default='verticalSensitivity.png',
                  help='Output PNG path (default: verticalSensitivity.png)')
    p.add_argument('--show', action='store_true', help='Display interactively')
    args = p.parse_args()

    gamma = np.linspace(0, 360, 721)

    fig, axes = plt.subplots(len(REGIONS), 2, figsize=(12, 11), sharex=True)

    for row, (label, relGeojsonA, relGeojsonD) in enumerate(REGIONS):
        geojsonA = os.path.join(args.projectDir, relGeojsonA)
        geojsonD = os.path.join(args.projectDir, relGeojsonD)
        H_A, psi_A, latA, lonA = headingAndIncidence(geojsonA)
        H_D, psi_D, latD, lonD = headingAndIncidence(geojsonD)
        alpha_deg = H_A - H_D
        H_mean = (H_A + H_D) / 2.0
        psi_mean = (psi_A + psi_D) / 2.0
        drToVelScale = 365.25 / (REPEAT_DAYS * np.sin(np.radians(psi_mean)))

        cx, cy = latLonToPS3413(0.5 * (latA + latD), 0.5 * (lonA + lonD))
        phi_real = np.degrees(np.arctan2(-cy, -cx)) % 360
        gamma_real = (phi_real - H_mean) % 360

        print(f'{label}: H_A={H_A:.1f} H_D={H_D:.1f} alpha={alpha_deg:.1f} '
              f'psi_mean={psi_mean:.1f} gamma_real={gamma_real:.1f} '
              f'drToVelScale={drToVelScale:.1f} (m/yr per m of {REPEAT_DAYS:.0f}-day dr error)')

        dvx = np.empty_like(gamma)
        dvy = np.empty_like(gamma)
        svx = np.empty_like(gamma)
        svy = np.empty_like(gamma)
        for i, g in enumerate(gamma):
            dvx[i], dvy[i], rnx, rny = sensitivities(alpha_deg, g, psi_A, psi_D)
            svx[i] = rnx * drToVelScale
            svy[i] = rny * drToVelScale

        axL, axR = axes[row, 0], axes[row, 1]

        # Left: independent dr-error sensitivity (primary interest)
        axL.plot(gamma, svx, label=r'$\sigma_{vx}$', color='C0')
        axL.plot(gamma, svy, label=r'$\sigma_{vy}$', color='C1')
        axL.set_ylabel(f'{label}\n(m/yr per m of {REPEAT_DAYS:.0f}-day $dr$ error)')
        axL.grid(alpha=0.3)
        axL.axvline(gamma_real, color='k', lw=1.2, ls='--',
                   label=f'real crossing ($\\gamma$={gamma_real:.0f}°)')

        # Right: common-mode vz sensitivity
        axR.plot(gamma, dvx, label=r'$\Delta v_x$', color='C0')
        axR.plot(gamma, dvy, label=r'$\Delta v_y$', color='C1')
        axR.set_ylabel('(per m/yr $v_z^{SMB}$)')
        axR.grid(alpha=0.3)
        axR.axvline(gamma_real, color='k', lw=1.2, ls='--')

        title = (f'asc H={H_A:.1f}°  x  desc H={H_D:.1f}°   '
                f'$\\alpha$={alpha_deg:.1f}°  $\\psi$={psi_mean:.1f}°')
        axL.set_title(title, fontsize=9)
        axR.set_title(title, fontsize=9)
        if row == 0:
            axL.legend(fontsize=8, loc='upper right')
            axR.legend(fontsize=8, loc='upper right')

        # Secondary top x-axis: actual longitude. lon = gamma + H_mean - 135
        # (verified numerically against the real EPSG:3413 transform -- 135 =
        # 90 + 45 comes from EPSG:3413's Longitude-of-origin=-45 baked into a
        # fixed linear offset). Deliberately UNWRAPPED (no mod-360) here --
        # secondary_xaxis needs a monotonic transform to place ticks
        # correctly; wrapping at +/-180 makes gammaToLon(0) == gammaToLon(360)
        # and confuses matplotlib's tick placement (collapses the axis to a
        # near-degenerate range). The wrap is applied only in the tick label
        # formatter below, for conventional +/-180 display.
        def gammaToLon(g, _Hm=H_mean):
            return g + _Hm - 135.0

        def lonToGamma(lon, _Hm=H_mean):
            return lon - _Hm + 135.0

        def wrappedLonFormatter(val, pos):
            return f'{((val + 180.0) % 360.0) - 180.0:.0f}'

        secL = axL.secondary_xaxis('top', functions=(gammaToLon, lonToGamma))
        secR = axR.secondary_xaxis('top', functions=(gammaToLon, lonToGamma))
        secL.xaxis.set_major_formatter(wrappedLonFormatter)
        secR.xaxis.set_major_formatter(wrappedLonFormatter)
        if row == 0:
            secL.set_xlabel('Longitude (deg)', fontsize=8)
            secR.set_xlabel('Longitude (deg)', fontsize=8)
        secL.tick_params(labelsize=7)
        secR.tick_params(labelsize=7)

    axes[0, 0].set_title('Independent dr-error sensitivity\n' + axes[0, 0].get_title(), fontsize=9)
    axes[0, 1].set_title('Uncompensated-$v_z$ sensitivity\n' + axes[0, 1].get_title(), fontsize=9)
    for ax in axes[-1, :]:
        ax.set_xlabel(r'$\gamma = \phi - H_{mean}$ (deg)')
    for ax in axes.flat:
        ax.set_xticks([0, 90, 180, 270, 360])

    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f'Written: {args.output}')
    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
