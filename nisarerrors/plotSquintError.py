#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plotSquintError.py

Quantifies the error in mosaic3d's computed (vx, vy) from ignoring residual
squint, at each of plotVerticalSensitivity's three real crossing-point
locations (northern/central/southern Greenland).

REVISED (2026-06-26): the original version of this script modeled squint via
a generalized nonzero-Doppler Newton solve (a *time*-shift at the pixel's
azimuth position). That model is WRONG -- confirmed directly against real
NISAR processor output. The image's azimuth-time assignment genuinely IS at
zero-Doppler (the time-solve in common/llToImageNew.c is correct and
untouched); residual squint instead shows up as a direct ANGULAR offset of
the true line-of-sight from the idealized zero-squint cross-track direction,
AT THAT SAME (unshifted) time/position. Confirmed via
science/LSAR/RUNW/metadata/geolocationGrid's losUnitVectorX/Y and
alongTrackUnitVectorX/Y cubes in real RUNW HDF5 files: the angle between them
is consistently ~91.5-91.9 deg, not exactly 90 -- and nisarBaseHDF.py's own
computeAngles() docstring confirms "the dataCube elevation angle includes the
squint" (i.e. the official processor's geolocation grid already reflects the
true squinted geometry, unlike GrIMP's own zero-Doppler-plane calculation).

Measured real squint at each region's real RUNW products (scene-center,
mid-height-level), ascending vs descending:

    Region    squint_asc (deg)   squint_desc (deg)
    North     1.658              1.494
    Central   1.667              1.481
    South     1.674              1.482

Ascending is consistently ~0.18 deg higher than descending across all three
regions -- a real, systematic asymmetry (not noise), plausibly related to
look/flight-direction-dependent yaw-steering-law residuals. Within one frame,
squint varies nearly linearly with slant range (~1.5-1.9 deg across one
swath) and is nearly constant in azimuth and height (see
Documents/plotSquintError.md for the full characterization, including
frame-to-frame stability and the open question of longer-track,
higher-than-linear behavior for e.g. Antarctic tracks spanning a much wider
heading/latitude range than these Greenland examples).

Mechanism: mosaic3d builds its solving matrix from the assumed zero-squint
headings, A0 = computeA(H_A(0), H_D(0)). The true heading is
H_true = H_0 + squint (a direct additive correction, NOT a time-shift),
giving A_true = computeA(H_A(0)+squint_A, H_D(0)+squint_D). For a true
velocity v_true, mosaic3d computes
v_computed = A0 @ A_true^-1 @ v_true = M @ v_true.

Two figures are generated:
  1. (--output) %error in vx/vy at theta=0/90/45 deg (see below), swept over
     mean squint angle (deg, x-axis) -- squint_A and squint_D are scaled
     together so each region's own real measured ratio between them is
     preserved at every point on the curve, up to 1.5x the real values.
  2. (--outputByDirection) at the real measured squint, the
     error decomposed as a function of true flow ORIENTATION theta (0-360
     deg). Left column: speed % error, plus vx/vy % error (each component's
     (computed - true) error, normalized by the true speed -- same
     normalization as speed % error, well-defined even where a component
     passes through zero) -- all three oscillate with theta (period 180
     deg). Right column: direction error in degrees (nearly CONSTANT
     regardless of theta, ~1.5-1.7 deg for all 3 regions). The near-constancy
     of the direction error is because M = A0 @ A_true^-1 is dominated by its
     antisymmetric part (M[0,1] ~= -M[1,0], much larger than M[0,0]-1 or
     M[1,1]-1) -- i.e. M is close to a pure small-angle rotation matrix, and a
     rotation rotates every input vector by the same angle regardless of its
     direction. See Documents/plotSquintError.md for the full writeup.

Usage:
    plotSquintError [-o OUTPUT.png] [--outputByDirection OUTPUT2.png] [--show]
"""
import argparse
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from nisarerrors.plotVerticalSensitivity import (
    computeA, headingAndIncidence, latLonToPS3413, REGIONS, DEFAULT_PROJECT_DIR)

SCALE_RANGE = (0.0, 1.5)
N_SCALE_POINTS = 151
N_THETA_POINTS = 721

# Real RUNW HDF5 paths for each region's ascending/descending sub-frame used
# to measure real squint (representative single sub-frame per region/image;
# see module docstring for the measured values and Documents/plotSquintError.md
# for the within-frame and frame-to-frame characterization).
SQUINT_SOURCE_H5 = {
    'Northern Greenland': (
        'track-58/source/NISAR_L1_PR_RUNW_015_058_A_037_016_7700_SH_20260314T051606_20260314T051641_20260326T051606_20260326T051641_P05012_N_P_J_001.h5',
        'track-97/source/NISAR_L1_PR_RUNW_015_097_D_053_016_7700_SH_20260316T222033_20260316T222110_20260328T222034_20260328T222110_P05012_N_P_J_001.h5'),
    'Central Greenland': (
        'track-131/source/NISAR_L1_PR_RUNW_007_131_A_037_008_7700_SH_20251213T064734_20251213T064815_20251225T064734_20251225T064811_X05010_N_P_J_001.h5',
        'track-11/source/NISAR_L1_PR_RUNW_004_011_D_055_006_7700_SH_20251029T231136_20251029T231213_20251122T231136_20251122T231213_X05010_N_P_J_001.h5'),
    'Southern Greenland': (
        'track-45/source/NISAR_L1_PR_RUNW_015_045_A_034_016_7700_SH_20260313T073547_20260313T073619_20260325T073547_20260325T073619_P05012_N_P_J_001.h5',
        'track-126/source/NISAR_L1_PR_RUNW_015_126_D_055_016_7700_SH_20260318T223823_20260318T223900_20260330T223823_20260330T223901_P05012_N_P_J_001.h5'),
}


def measureSquint(h5path):
    """Real squint (deg) at scene center (mid height/azimuth/range), measured
    as the deviation of the angle between losUnitVector and
    alongTrackUnitVector from the idealized 90 deg (broadside). Both vectors
    come directly from the official processor's geolocationGrid cube -- no
    geometric modeling/assumptions, see module docstring."""
    h5 = h5py.File(h5path, 'r')
    grp = h5['science']['LSAR']['RUNW']['metadata']['geolocationGrid']
    losX, losY = grp['losUnitVectorX'][:], grp['losUnitVectorY'][:]
    atX, atY = grp['alongTrackUnitVectorX'][:], grp['alongTrackUnitVectorY'][:]
    losHeading = np.degrees(np.arctan2(losX, losY)) % 360
    atHeading = np.degrees(np.arctan2(atX, atY)) % 360
    diff = ((losHeading - atHeading + 180) % 360) - 180
    squint = diff - 90.0
    nh, na, nr = squint.shape
    h5.close()
    return float(squint[nh // 2, na // 2, nr // 2])


def regionGeometry(label, relGeojsonA, relGeojsonD, projectDir):
    """Assumed-zero-squint matrix A0 and each region's real measured
    (squint_A, squint_D), plus the inputs needed to rebuild A_true at an
    arbitrary scale factor."""
    geojsonA = os.path.join(projectDir, relGeojsonA)
    geojsonD = os.path.join(projectDir, relGeojsonD)
    H_A0, psi_A, latA, lonA = headingAndIncidence(geojsonA)
    H_D0, psi_D, latD, lonD = headingAndIncidence(geojsonD)
    cx, cy = latLonToPS3413(0.5 * (latA + latD), 0.5 * (lonA + lonD))
    phi_real = np.degrees(np.arctan2(-cy, -cx)) % 360
    beta0 = np.radians(phi_real - H_A0)
    alpha0 = np.radians(H_A0 - H_D0)
    A0 = computeA(alpha0, beta0)

    relH5A, relH5D = SQUINT_SOURCE_H5[label]
    squintA_real = measureSquint(os.path.join(projectDir, relH5A))
    squintD_real = measureSquint(os.path.join(projectDir, relH5D))
    return A0, H_A0, H_D0, phi_real, squintA_real, squintD_real


def trueMatrix(H_A0, H_D0, phi_real, squintA, squintD, scale):
    H_A = H_A0 + scale * squintA
    H_D = H_D0 + scale * squintD
    beta_t = np.radians(phi_real - H_A)
    alpha_t = np.radians(H_A - H_D)
    return computeA(alpha_t, beta_t)


def errorAtTheta(M, thetaDeg):
    """Speed/vx/vy % error and direction error (deg) of M @ vhat(thetaDeg), for a
    single true-flow orientation. vx/vy % error is (computed - true) component,
    normalized by the true SPEED (always 1.0, since vhat is a unit vector) rather
    than by the true component itself -- the same normalization speed % error
    already uses, and one that stays well-defined as either true component passes
    through zero. At thetaDeg=0 this reduces to (M[0,0]-1)*100 (Figure 1's pure-x
    case); at thetaDeg=90, (M[1,1]-1)*100 (Figure 1's pure-y case)."""
    th = np.radians(thetaDeg)
    vhat = np.array([np.cos(th), np.sin(th)])
    vcomp = M @ vhat
    speedErrPct = (np.linalg.norm(vcomp) - 1.0) * 100.0
    vxErrPct = (vcomp[0] - vhat[0]) * 100.0
    vyErrPct = (vcomp[1] - vhat[1]) * 100.0
    diff = np.degrees(np.arctan2(vcomp[1], vcomp[0]) - th)
    angErrDeg = ((diff + 180) % 360) - 180   # wrap to (-180, 180]
    return speedErrPct, vxErrPct, vyErrPct, angErrDeg


def errorByDirection(M, nTheta=N_THETA_POINTS):
    """Speed/vx/vy % error and direction error (deg) of M @ vhat(theta), swept
    over every true-flow orientation theta (0-360 deg) -- see errorAtTheta."""
    thetaDeg = np.linspace(0, 360, nTheta)
    speedErrPct = np.empty_like(thetaDeg)
    vxErrPct = np.empty_like(thetaDeg)
    vyErrPct = np.empty_like(thetaDeg)
    angErrDeg = np.empty_like(thetaDeg)
    for i, th in enumerate(thetaDeg):
        speedErrPct[i], vxErrPct[i], vyErrPct[i], angErrDeg[i] = errorAtTheta(M, th)
    return thetaDeg, speedErrPct, vxErrPct, vyErrPct, angErrDeg


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--projectDir', default=DEFAULT_PROJECT_DIR,
                  help=f'Root directory containing the track-N/ example data '
                       f'(default: {DEFAULT_PROJECT_DIR})')
    p.add_argument('-o', '--output', default='squintError.png',
                  help='Output PNG path for the scale-sweep figure (default: squintError.png)')
    p.add_argument('--outputByDirection', default='squintErrorByDirection.png',
                  help='Output PNG path for the orientation-dependence figure '
                       '(default: squintErrorByDirection.png)')
    p.add_argument('--show', action='store_true', help='Display interactively')
    args = p.parse_args()

    scales = np.linspace(*SCALE_RANGE, N_SCALE_POINTS)

    fig1, axes1 = plt.subplots(len(REGIONS), 1, figsize=(8, 11), sharex=True)
    fig2, axes2 = plt.subplots(len(REGIONS), 2, figsize=(11, 9))

    for row, (label, relGeojsonA, relGeojsonD) in enumerate(REGIONS):
        A0, H_A0, H_D0, phi_real, squintA_real, squintD_real = regionGeometry(
            label, relGeojsonA, relGeojsonD, args.projectDir)

        # Sweep both squint_A and squint_D together, scaled by the same factor,
        # so their real measured ratio (the actual asc/desc asymmetry) is
        # preserved at every point -- but plot directly against the resulting
        # mean squint angle in degrees rather than the dimensionless factor
        # itself, since that's what's physically varying and needs no
        # secondary axis or title lookup to interpret.
        meanSquintReal = 0.5 * (squintA_real + squintD_real)
        meanSquintDeg = scales * meanSquintReal

        pctErrVx = np.empty_like(scales)
        pctErrVy = np.empty_like(scales)
        pctErrVx45 = np.empty_like(scales)
        pctErrVy45 = np.empty_like(scales)
        for i, scale in enumerate(scales):
            A_true = trueMatrix(H_A0, H_D0, phi_real, squintA_real, squintD_real, scale)
            M = A0 @ np.linalg.inv(A_true)
            _, pctErrVx[i], _, _ = errorAtTheta(M, 0.0)
            _, _, pctErrVy[i], _ = errorAtTheta(M, 90.0)
            _, pctErrVx45[i], pctErrVy45[i], _ = errorAtTheta(M, 45.0)

        print(f'{label}: squint_asc={squintA_real:.4f} deg  squint_desc={squintD_real:.4f} deg  '
              f'%error_vx range=[{pctErrVx.min():.3f},{pctErrVx.max():.3f}]  '
              f'%error_vy range=[{pctErrVy.min():.3f},{pctErrVy.max():.3f}]  '
              f'%error_vx(45deg) range=[{pctErrVx45.min():.3f},{pctErrVx45.max():.3f}]  '
              f'%error_vy(45deg) range=[{pctErrVy45.min():.3f},{pctErrVy45.max():.3f}]')

        ax = axes1[row]
        ax.plot(meanSquintDeg, pctErrVx, label=r'$\%$error $v_x$ ($\theta{=}0°$)', color='C0')
        ax.plot(meanSquintDeg, pctErrVy, label=r'$\%$error $v_y$ ($\theta{=}90°$)', color='C1')
        ax.plot(meanSquintDeg, pctErrVx45, label=r'$\%$error $v_x$ ($\theta{=}45°$)', color='C0', ls='--')
        ax.plot(meanSquintDeg, pctErrVy45, label=r'$\%$error $v_y$ ($\theta{=}45°$)', color='C1', ls='--')
        ax.axhline(0.0, color='gray', lw=0.5, ls=':')
        ax.axvline(meanSquintReal, color='k', lw=1.2, ls='--',
                  label='real measured squint' if row == 0 else None)
        ax.set_ylabel(f'{label}\n% error')
        ax.grid(alpha=0.3)
        if row == 0:
            ax.legend(fontsize=8, loc='upper left')
        ax.set_title(f'real squint: asc={squintA_real:.2f}°  desc={squintD_real:.2f}°  '
                    f'(mean={meanSquintReal:.2f}°, dashed line)', fontsize=9)

        # Second figure: error decomposed by true-flow orientation, at the real measured squint
        M_real = A0 @ np.linalg.inv(
            trueMatrix(H_A0, H_D0, phi_real, squintA_real, squintD_real, 1.0))
        thetaDeg, speedErrPct, vxErrPct, vyErrPct, angErrDeg = errorByDirection(M_real)

        print(f'  speed %err=[{speedErrPct.min():.3f},{speedErrPct.max():.3f}]  '
              f'vx %err=[{vxErrPct.min():.3f},{vxErrPct.max():.3f}]  '
              f'vy %err=[{vyErrPct.min():.3f},{vyErrPct.max():.3f}]  '
              f'direction err=[{angErrDeg.min():.3f},{angErrDeg.max():.3f}] deg')

        ax1, ax2 = axes2[row]
        ax1.plot(thetaDeg, speedErrPct, label='speed' if row == 0 else None, color='C0')
        ax1.plot(thetaDeg, vxErrPct, label=r'$v_x$' if row == 0 else None, color='C1')
        ax1.plot(thetaDeg, vyErrPct, label=r'$v_y$' if row == 0 else None, color='C3')
        ax1.axhline(0, color='gray', lw=0.5, ls=':')
        ax1.set_ylabel(f'{label}\n% error')
        ax1.grid(alpha=0.3)
        if row == 0:
            ax1.legend(fontsize=8, loc='upper right')
        ax2.plot(thetaDeg, angErrDeg, color='C2')
        ax2.axhline(0, color='gray', lw=0.5, ls=':')
        ax2.set_ylabel('direction error (deg)')
        ax2.grid(alpha=0.3)

    axes1[-1].set_xlabel('Mean squint angle, (squint_A+squint_D)/2 (deg) -- '
                         'squint_A/squint_D scaled together, real ratio preserved')
    fig1.tight_layout()
    fig1.savefig(args.output, dpi=150)
    print(f'Written: {args.output}')

    axes2[-1, 0].set_xlabel('True flow direction theta (deg, ccw from x-axis)')
    axes2[-1, 1].set_xlabel('True flow direction theta (deg, ccw from x-axis)')
    fig2.suptitle('Velocity error vs true flow direction, at real measured squint')
    fig2.tight_layout()
    fig2.savefig(args.outputByDirection, dpi=150)
    print(f'Written: {args.outputByDirection}')

    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
