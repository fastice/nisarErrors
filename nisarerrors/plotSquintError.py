#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plotSquintError.py

mosaic3d always assumes zero-Doppler (broadside) geometry when computing
heading (common/computeHeading.c), via a Newton solve that targets exactly
zero Doppler (utilities.geodatrxa.geodatrxa.llzPtToRA()). Real acquisitions
can have a small residual squint (nonzero Doppler) even after zero-Doppler
*processing* -- this script quantifies the resulting error in mosaic3d's
computed (vx, vy) at each of plotVerticalSensitivity.py's three real
crossing-point locations (the same REGIONS, evaluated only at each region's
already-marked real gamma -- not swept over pixel position).

Mechanism: mosaic3d builds its solving matrix from the ASSUMED zero-squint
headings, A0 = computeA(H_A(0), H_D(0)). The TRUE measurement reflects
whatever the ACTUAL (possibly squinted) heading was,
A_true = computeA(H_A(f_dopA), H_D(f_dopD)). For a true velocity v_true,
mosaic3d computes v_computed = A0 @ A_true^-1 @ v_true = M @ v_true. Since M
is linear, the resulting fractional error is independent of |v_true| -- only
its direction matters in general. We report each component's error using ITS
OWN natural reference case (true velocity purely along that axis), which is
well-defined with no extra direction assumption and no division-by-zero:

    %error_vx = (M[0,0] - 1) * 100   (true velocity = pure x-hat)
    %error_vy = (M[1,1] - 1) * 100   (true velocity = pure y-hat)

Heading at a nonzero target Doppler is computed RIGOROUSLY (not via a small-
angle approximation -- confirmed analytically inadequate, off by a factor of
~4 numerically) by generalizing geodatrxa.llzPtToRA()'s Newton solve to a
nonzero target:

    zero-Doppler (llzPtToRA):       (T-S).V = 0
    target f_dop (this script):     (T-S).V = f_dop*lambda*R/2

derived from the standard f_dop = -(2/lambda)*dR/dt = (2/lambda)*(T-S).V/R
relation. Reuses only geodatrxa's Doppler-agnostic primitives (lltoecef,
interpPos, interpVel, rNearSLP, slpRg, prf) -- geodatrxa.py itself is not
modified.

The squint is applied as the SAME squint/PRF ratio to both the ascending and
descending image (each using its own real PRF) -- i.e. a shared fractional
residual squint, not independently-varying squints per image.

Note: this captures only the heading (A-matrix) error -- a separate, smaller
effect (~5-6x less, confirmed numerically) on the psi-dependent phase/offset-
to-velocity-units scaling is NOT included here, since computeA() itself does
not depend on psi at all (verified against common/initRoutines.c) and whether
mosaic3d's geodat-metadata psi already reflects the true acquisition geometry
is ambiguous without knowing the NISAR processor's internals.

Usage:
    plotSquintError [--projectDir DIR] [-o OUTPUT.png] [--show]
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from utilities import geodatrxa
from nisarerrors.plotVerticalSensitivity import (
    computeA, headingAndIncidence, latLonToPS3413, REGIONS, DEFAULT_PROJECT_DIR)

SQUINT_RATIO_RANGE = (-2.0, 2.0)
N_SQUINT_POINTS = 161


def headingAtDoppler(g, lat, lon, f_dop):
    """Cross-track heading (deg from north, CW) at (lat, lon) for an already-
    loaded geodatrxa object g, targeting Doppler frequency f_dop (Hz) instead
    of geodatrxa.llzPtToRA()'s hardcoded zero. Same structure/formula as
    headingAndIncidence() in plotVerticalSensitivity.py, but with the
    generalized Newton solve below in place of g.llzPtToRA().
    """
    def squintedRA(latPt, lonPt, z, initT=None):
        tPt = np.array(g.lltoecef(latPt, lonPt, z))
        if initT is None:
            initT = g.t0 + 0.5 * g.na * g.nla / g.prf
        myTime = initT
        for _ in range(50):
            sPt = np.array(g.interpPos(myTime))
            vPt = np.array(g.interpVel(myTime))
            dr = tPt - sPt
            R = np.sqrt(np.dot(dr, dr))
            target = f_dop * g.wavelength * R / 2.0
            df = np.dot(dr, vPt) - target
            C1 = -np.dot(vPt, vPt)
            myTime -= df / C1
            if np.abs(df / C1) < 1e-5:
                break
        sPt = np.array(g.interpPos(myTime))
        dr = tPt - sPt
        r = (np.sqrt(np.dot(dr, dr)) - g.rNearSLP) / g.slpRg
        az = (myTime - g.t0) * g.prf
        return r, az, myTime

    dlat = -0.0025 if g.lookdir == 'left' else 0.0025
    r1, a1, t1 = squintedRA(lat - dlat, lon, 0.0)
    r2, a2, t2 = squintedRA(lat + dlat, lon, 0.0, initT=t1)
    da = (a2 - a1) * g.slpAz
    Re = g.earthRadLatm(np.radians(lat))
    ReH1, ReH2 = g.ReH(t1), g.ReH(t2)
    R1 = r1 * g.slpRg + g.rNearSLP
    R2 = r2 * g.slpRg + g.rNearSLP
    rho1 = np.arccos((R1 ** 2 - ReH1 ** 2 - Re ** 2) / (-2.0 * Re * ReH1))
    rho2 = np.arccos((R2 ** 2 - ReH2 ** 2 - Re ** 2) / (-2.0 * Re * ReH2))
    dgr = Re * rho2 - Re * rho1
    hAngle = np.arctan2(da, dgr) if g.lookdir == 'right' else np.arctan2(da, -dgr)
    return np.degrees(hAngle) % 360


def squintAngleDeg(ratio, prf, wavelength, vSat):
    """True antenna squint angle (deg) for a given squint/PRF ratio -- the
    physical angle between the LOS and broadside, NOT the resulting heading
    shift (which is a separate, smaller, non-proportional effect -- see
    module docstring)."""
    fdop = ratio * prf
    return np.degrees(np.arcsin(np.clip(fdop * wavelength / (2.0 * vSat), -1, 1)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--projectDir', default=DEFAULT_PROJECT_DIR,
                  help=f'Root directory containing the track-N/ example geodats '
                       f'(default: {DEFAULT_PROJECT_DIR})')
    p.add_argument('-o', '--output', default='squintError.png',
                  help='Output PNG path (default: squintError.png)')
    p.add_argument('--show', action='store_true', help='Display interactively')
    args = p.parse_args()

    ratios = np.linspace(*SQUINT_RATIO_RANGE, N_SQUINT_POINTS)

    fig, axes = plt.subplots(len(REGIONS), 1, figsize=(8, 11), sharex=True)

    for row, (label, relGeojsonA, relGeojsonD) in enumerate(REGIONS):
        geojsonA = os.path.join(args.projectDir, relGeojsonA)
        geojsonD = os.path.join(args.projectDir, relGeojsonD)

        # Zero-squint baseline (same as plotVerticalSensitivity.py)
        H_A0, psi_A, latA, lonA = headingAndIncidence(geojsonA)
        H_D0, psi_D, latD, lonD = headingAndIncidence(geojsonD)
        cx, cy = latLonToPS3413(0.5 * (latA + latD), 0.5 * (lonA + lonD))
        phi_real = np.degrees(np.arctan2(-cy, -cx)) % 360
        beta0 = np.radians(phi_real - H_A0)
        alpha0 = np.radians(H_A0 - H_D0)
        A0 = computeA(alpha0, beta0)

        gA = geodatrxa(); gA.readGeojson(geojsonA)
        gD = geodatrxa(); gD.readGeojson(geojsonD)
        vSatA = np.linalg.norm(gA.velocity[len(gA.velocity) // 2])

        pctErrVx = np.empty_like(ratios)
        pctErrVy = np.empty_like(ratios)
        for i, ratio in enumerate(ratios):
            H_A = headingAtDoppler(gA, latA, lonA, ratio * gA.prf)
            H_D = headingAtDoppler(gD, latD, lonD, ratio * gD.prf)
            beta_t = np.radians(phi_real - H_A)
            alpha_t = np.radians(H_A - H_D)
            A_true = computeA(alpha_t, beta_t)
            M = A0 @ np.linalg.inv(A_true)
            pctErrVx[i] = (M[0, 0] - 1.0) * 100.0
            pctErrVy[i] = (M[1, 1] - 1.0) * 100.0

        print(f'{label}: PRF={gA.prf:.1f} Hz  lambda={gA.wavelength:.4f} m  '
              f'Vsat={vSatA:.1f} m/s  gamma_real={(phi_real - (H_A0+H_D0)/2.0)%360:.1f} '
              f'%error_vx range=[{pctErrVx.min():.2f},{pctErrVx.max():.2f}]  '
              f'%error_vy range=[{pctErrVy.min():.2f},{pctErrVy.max():.2f}]')

        ax = axes[row]
        ax.plot(ratios, pctErrVx, label=r'$\%$error $v_x$', color='C0')
        ax.plot(ratios, pctErrVy, label=r'$\%$error $v_y$', color='C1')
        ax.axhline(0.0, color='gray', lw=0.5, ls=':')
        ax.set_ylabel(f'{label}\n% error')
        ax.grid(alpha=0.3)
        if row == 0:
            ax.legend(fontsize=8, loc='upper right')
        ax.set_title(f'PRF={gA.prf:.0f} Hz  $\\lambda$={gA.wavelength:.3f} m  '
                    f'$V_{{sat}}$={vSatA:.0f} m/s', fontsize=9)

        # Secondary top x-axis: true antenna squint angle (deg), region-specific
        # (depends on this region's real PRF/wavelength/Vsat).
        def ratioToSquintDeg(r, _prf=gA.prf, _wl=gA.wavelength, _vs=vSatA):
            return squintAngleDeg(r, _prf, _wl, _vs)

        def squintDegToRatio(deg, _prf=gA.prf, _wl=gA.wavelength, _vs=vSatA):
            return 2.0 * _vs * np.sin(np.radians(deg)) / (_prf * _wl)

        sec = ax.secondary_xaxis('top', functions=(ratioToSquintDeg, squintDegToRatio))
        sec.set_xlabel('True antenna squint angle (deg)' if row == 0 else '', fontsize=8)
        sec.tick_params(labelsize=7)

    axes[-1].set_xlabel('squint / PRF')
    fig.tight_layout()
    fig.savefig(args.output, dpi=150)
    print(f'Written: {args.output}')
    if args.show:
        plt.show()


if __name__ == '__main__':
    main()
