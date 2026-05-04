#!/usr/bin/env python
"""
Validate mixed wide+delensing LAT noise maps produced by wide-delens_mix.py.

For a given DELENS_FRACTION f the expected per-pixel noise variance is:

  Wide-only pixels   : var_comb(p) = var_wide(p)   / (1-f)
  Delens-only pixels : var_comb(p) = var_delens(p) /  f
  Both-survey pixels : 1/var_comb(p) = (1-f)/var_wide(p) + f/var_delens(p)

This script:
  1. Builds the expected combined variance map from the stored variance maps.
  2. Plots a mollview of the first noise realization for each channel.
  3. Estimates the average power spectra (hits-weighted anafast + NMT-decoupled).
  4. Also computes a homogenized NMT spectrum restricted to the analysis mask.
  5. Compares measured white-noise level and 1/f shape against theory.
  6. Reports measured vs expected depths.
"""

import os
import yaml
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import cmocean
import pymaster as nmt
from tqdm import tqdm

# ── Configuration ────────────────────────────────────────────────────────────

NSIMS_FOR_AVG = 5

# Must match the value used when running wide-delens_mix.py
DELENS_FRACTION = 0.5

WIDE_YAML   = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/instr_params_lat-wide_channels.yaml'
DELENS_YAML = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/instr_params_lat-delens_channels.yaml'

BASE_OUTPUT_DIR    = '/pscratch/sd/s/shamikg/so_mapbased_noise/output'
VARIANCE_MAP_DIR   = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/variance_maps'
ANALYSIS_MASK_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/so_lat_deep_wide_analysis_2048.fits'

LMAX = 4500

# First-5-char channel prefix → frequency code in variance-map filenames
CHANNEL_FREQ_MAP = {
    'LF027': 'f027', 'LF039': 'f039',
    'MF093': 'f093', 'MF145': 'f145',
    'HF225': 'f225', 'HF280': 'f280',
}


# ── Helpers ──────────────────────────────────────────────────────────────────

def load_variance_map(ch_name, tel_yrs, survey, variance_map_dir):
    """Load the per-pixel noise variance map for a channel / survey / tel-yrs."""
    freq_code = CHANNEL_FREQ_MAP[ch_name[:5]]
    if survey == 'wide':
        fname = f"so_lat_wide_1_el_{freq_code}_noise_var_{tel_yrs:.2f}tel-yrs.fits"
    elif survey == 'delens_wide':
        fname = f"so_lat_delens_wide_{freq_code}_noise_var_{tel_yrs:.2f}tel-yrs.fits"
    else:
        raise ValueError(f"Unknown survey: {survey!r}")
    fpath = os.path.join(variance_map_dir, fname)
    if not os.path.exists(fpath):
        raise FileNotFoundError(f"Variance map not found: {fpath}")
    return hp.read_map(fpath, verbose=False)


def combined_variance_map(var_wide, var_delens, f):
    """
    Compute the per-pixel noise variance of the ILC-combined map for split fraction f.

    Derivation:
      After rescaling each survey's budget, the new variances are
        var_w_new = var_wide   / (1-f)
        var_d_new = var_delens / f
      ILC combination of two independent streams minimises variance:
        1/var_comb = 1/var_w_new + 1/var_d_new
                   = (1-f)/var_wide + f/var_delens
    """
    a = 1.0 - f   # wide fraction
    b = f          # delensing fraction

    w_in = var_wide   > 0
    d_in = var_delens > 0
    var_comb = np.zeros(len(var_wide))

    # Wide-only pixels
    w_only = w_in & ~d_in
    if np.any(w_only) and a > 0:
        var_comb[w_only] = var_wide[w_only] / a

    # Delensing-only pixels
    d_only = ~w_in & d_in
    if np.any(d_only) and b > 0:
        var_comb[d_only] = var_delens[d_only] / b

    # Pixels in both footprints — ILC
    both = w_in & d_in
    if np.any(both):
        if a == 0:
            var_comb[both] = var_delens[both]
        elif b == 0:
            var_comb[both] = var_wide[both]
        else:
            vw = var_wide[both]
            vd = var_delens[both]
            var_comb[both] = vw * vd / (a * vd + b * vw)

    return var_comb


def Nl2uKarcmin(Nl):
    """Convert noise power Nl [uK^2 sr] to noise depth [uK-arcmin]."""
    return np.sqrt(Nl) / np.deg2rad(1. / 60.)


def get_theory_Nl_ip(Nl_white_I, Nl_white_P, params, lmax):
    """
    Build theory noise power spectrum curves including the 1/f component.

    N_l^I(l) = Nl_white_I * (1 + (l/l_knee_i)^alpha_i)
    N_l^P(l) = Nl_white_P * (1 + (l/l_knee_p)^alpha_p)
    """
    ells = np.arange(lmax + 1)
    factor_i = np.ones(lmax + 1)
    factor_i[2:] = 1. + (ells[2:] / params['ell_knee_i']) ** params['alpha_knee_i']
    factor_p = np.ones(lmax + 1)
    factor_p[2:] = 1. + (ells[2:] / params['ell_knee_p']) ** params['alpha_knee_p']
    return ells, Nl_white_I * factor_i, Nl_white_P * factor_p


def noise_filepath(noise_dir, ch_name, f, sim_idx, nside):
    """Return path to a mixed noise-map FITS file."""
    return os.path.join(
        noise_dir,
        f"sobs_lat-noise-mix_{ch_name}_f{f:.2f}_mc{sim_idx:03d}_nside{nside:04d}.fits",
    )


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    f = float(DELENS_FRACTION)
    assert 0.0 <= f <= 1.0, "DELENS_FRACTION must be in [0, 1]"

    with open(WIDE_YAML)   as fh: wide_params  = yaml.safe_load(fh)
    with open(DELENS_YAML) as fh: delens_params = yaml.safe_load(fh)

    # Load analysis mask — applied to restrict the 3rd power spectrum track
    print(f"Loading analysis mask: {ANALYSIS_MASK_FILE}")
    analysis_mask = hp.read_map(ANALYSIS_MASK_FILE, verbose=False)

    channels = list(wide_params.keys())
    if set(channels) != set(delens_params.keys()):
        print("Warning: channel lists differ between YAMLs — using intersection.")
        channels = sorted(set(channels) & set(delens_params.keys()))

    mix_noise_dir     = os.path.join(BASE_OUTPUT_DIR, f'lat-mix_f{f:.2f}')
    validation_folder = os.path.join(mix_noise_dir, 'validation')
    os.makedirs(validation_folder, exist_ok=True)

    print(f"Delensing fraction  f = {f:.2f}  "
          f"(wide: {1-f:.2f}, delensing: {f:.2f})")
    print(f"Reading noise maps from : {mix_noise_dir}")
    print(f"Saving validation to    : {validation_folder}")
    print(f"Channels ({len(channels)}): {channels}\n")

    for ch_name in channels:
        print(f"\n{'='*60}")
        print(f"Channel: {ch_name}")
        print(f"{'='*60}")

        wp    = wide_params[ch_name]
        dp    = delens_params[ch_name]
        nside = wp['nside']

        # ── Build combined variance map ───────────────────────────────────────
        var_wide   = load_variance_map(ch_name, wp['tel_yrs'], 'wide',        VARIANCE_MAP_DIR)
        var_delens = load_variance_map(ch_name, dp['tel_yrs'], 'delens_wide', VARIANCE_MAP_DIR)
        var_comb   = combined_variance_map(var_wide, var_delens, f)

        obs          = var_comb > 0
        binary_mask  = obs.astype(float)
        fsky_binary  = float(np.mean(binary_mask))

        # Inverse-variance weight (analogous to relhits): normalised to [0, 1]
        weight = np.zeros_like(var_comb)
        weight[obs] = 1.0 / var_comb[obs]
        weight /= weight.max()
        # Effective fsky for pseudo-Cl: mean(w^2) / mean(w)  (standard approx)
        fsky = float(np.mean(weight ** 2) / np.mean(weight))

        # hp.mollview(weight, cmap='coolwarm')
        # plt.savefig(f'{validation_folder}/weights.png')
        # plt.close()
        # exit()
        mask_apo     = nmt.mask_apodization(binary_mask, 1.0, apotype='C2')
        mask_apo_anl = nmt.mask_apodization(binary_mask * analysis_mask, 1.0, apotype='C2')

        print(f"  Wide   tel-yrs: {wp['tel_yrs']}  |  Delens tel-yrs: {dp['tel_yrs']}")
        print(f"  f_sky (binary): {fsky_binary:.4f}  |  f_sky (eff, inv-var wtd): {fsky:.4f}")

        # Expected white noise from variance maps (var_comb is in uK²-pixel units)
        # Conversion: depth [uK-arcmin] = 60 * sqrt(var_pix * pixarea_deg)
        # In _get_noise_variance_map: I is divided by sqrt(2), so depth_T = depth_P / sqrt(2)
        pixarea_deg = hp.nside2pixarea(nside, degrees=True)
        depth_exp_P = float(60. * np.sqrt(np.median(var_comb[weight == 1.]) * pixarea_deg))
        depth_exp_T = depth_exp_P / np.sqrt(2.)
        print(f"  Expected depth from var map  — T: {depth_exp_T:.2f}  P: {depth_exp_P:.2f}  uK-arcmin")

        # ── Mollview of first realization ─────────────────────────────────────
        noise_file_0 = noise_filepath(mix_noise_dir, ch_name, f, 0, nside)
        if not os.path.exists(noise_file_0):
            print(f"  Noise map not found: {noise_file_0} — skipping channel")
            continue

        noise_map_0 = hp.read_map(noise_file_0, field=None, verbose=False)

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        for i, (ax, title) in enumerate(zip(axes, ['I (Temperature)', 'Q', 'U'])):
            plt.sca(ax)
            valid = noise_map_0[i][obs]
            vmax  = float(np.std(valid)) * 3
            hp.mollview(
                noise_map_0[i],
                title=f'{ch_name}  f={f:.2f} — {title}',
                cmap=cmocean.cm.balance,
                min=-vmax, max=vmax,
                unit=r'$\mu K$',
                hold=True,
            )
            hp.graticule(lw=0.5, color='gray', ls='--', alpha=0.5)
        plt.savefig(
            os.path.join(validation_folder, f'mollview_{ch_name}_f{f:.2f}_mc000.png'),
            dpi=150, bbox_inches='tight',
        )
        plt.close()

        # ── Average power spectra ─────────────────────────────────────────────
        print(f"  Computing average power spectra over {NSIMS_FOR_AVG} sims...")

        bins    = nmt.NmtBin.from_lmax_linear(LMAX, 30)
        leff    = bins.get_effective_ells()
        nmt_wsp = nmt.NmtWorkspace()

        cl_sum_TT = np.zeros(LMAX + 1)
        cl_sum_EE = np.zeros(LMAX + 1)
        cl_sum_BB = np.zeros(LMAX + 1)
        cl_h_sum_TT = np.zeros(len(leff))
        cl_h_sum_EE = np.zeros(len(leff))
        cl_h_sum_BB = np.zeros(len(leff))
        # Analysis-mask-restricted homogenized spectra
        cl_ha_sum_TT = np.zeros(len(leff))
        cl_ha_sum_EE = np.zeros(len(leff))
        cl_ha_sum_BB = np.zeros(len(leff))
        nmt_wsp_anl  = nmt.NmtWorkspace()
        n_done = 0

        for sim_idx in tqdm(range(NSIMS_FOR_AVG), desc=f"  {ch_name}", ncols=120):
            noise_file = noise_filepath(mix_noise_dir, ch_name, f, sim_idx, nside)
            if not os.path.exists(noise_file):
                print(f"  Missing sim {sim_idx}, skipping...")
                continue

            nm = hp.read_map(noise_file, field=None, verbose=False)

            # Inv-variance weighted pseudo-Cl (quick estimate)
            cls = hp.anafast(nm * weight, lmax=LMAX, pol=True, nspec=3) / fsky
            cl_sum_TT += cls[0]
            cl_sum_EE += cls[1]
            cl_sum_BB += cls[2]

            # Homogenized map for NMT (multiply by sqrt(inv-var weight))
            # This makes the effective per-pixel variance ~constant before masking
            homo_map = nm * np.sqrt(weight)

            if sim_idx == 0:
                field_T = nmt.NmtField(
                    mask_apo, [homo_map[0]],
                    purify_e=False, purify_b=False,
                    lmax=LMAX, lmax_mask=LMAX,
                )
                nmt_wsp.compute_coupling_matrix(field_T, field_T, bins)
                field_T_anl = nmt.NmtField(
                    mask_apo_anl, [homo_map[0]],
                    purify_e=False, purify_b=False,
                    lmax=LMAX, lmax_mask=LMAX,
                )
                nmt_wsp_anl.compute_coupling_matrix(field_T_anl, field_T_anl, bins)

            cls_h = hp.anafast(homo_map * mask_apo, lmax=LMAX, pol=True, nspec=3)
            for arr, cl_row in zip(
                [cl_h_sum_TT, cl_h_sum_EE, cl_h_sum_BB], cls_h[:3]
            ):
                arr += nmt_wsp.decouple_cell(cl_row.reshape(1, -1))[0]

            # Analysis-mask-restricted homogenized spectra
            cls_ha = hp.anafast(homo_map * mask_apo_anl, lmax=LMAX, pol=True, nspec=3)
            for arr, cl_row in zip(
                [cl_ha_sum_TT, cl_ha_sum_EE, cl_ha_sum_BB], cls_ha[:3]
            ):
                arr += nmt_wsp_anl.decouple_cell(cl_row.reshape(1, -1))[0]

            n_done += 1

        if n_done == 0:
            print(f"  No simulations found for {ch_name}, skipping.")
            continue

        cl_avg_TT    = cl_sum_TT    / n_done
        cl_avg_EE    = cl_sum_EE    / n_done
        cl_avg_BB    = cl_sum_BB    / n_done
        cl_h_avg_TT  = cl_h_sum_TT  / n_done
        cl_h_avg_EE  = cl_h_sum_EE  / n_done
        cl_h_avg_BB  = cl_h_sum_BB  / n_done
        cl_ha_avg_TT = cl_ha_sum_TT / n_done
        cl_ha_avg_EE = cl_ha_sum_EE / n_done
        cl_ha_avg_BB = cl_ha_sum_BB / n_done

        # ── White-noise estimate from high-ell ────────────────────────────────
        idx_lo, idx_hi = 3000, min(4500, LMAX)
        Nl_white_TT_meas = float(np.mean(cl_avg_TT[idx_lo:idx_hi]))
        Nl_white_P_meas  = float(np.mean(
            [np.mean(cl_avg_EE[idx_lo:idx_hi]),
             np.mean(cl_avg_BB[idx_lo:idx_hi])]
        ))
        uKarcmin_TT = Nl2uKarcmin(Nl_white_TT_meas)
        uKarcmin_P  = Nl2uKarcmin(Nl_white_P_meas)

        print(f"  Measured  depth (ell {idx_lo}–{idx_hi}) — "
              f"T: {uKarcmin_TT:.2f}  P: {uKarcmin_P:.2f}  uK-arcmin")
        print(f"  P/T ratio: {uKarcmin_P/uKarcmin_TT:.3f}  (expected ~1.414)")

        # Theory curves from fitted white noise floor
        ells, Nl_theory_T, Nl_theory_P = get_theory_Nl_ip(
            Nl_white_TT_meas, Nl_white_P_meas, wp, LMAX,
        )

        # ── Plot ──────────────────────────────────────────────────────────────
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        specs = [
            ('TT', cl_avg_TT, cl_h_avg_TT, cl_ha_avg_TT, Nl_theory_T, 'tab:blue'),
            ('EE', cl_avg_EE, cl_h_avg_EE, cl_ha_avg_EE, Nl_theory_P, 'tab:orange'),
            ('BB', cl_avg_BB, cl_h_avg_BB, cl_ha_avg_BB, Nl_theory_P, 'tab:green'),
        ]
        for ax, (label, cl_meas, cl_h_meas, cl_ha_meas, Nl_th, color) in zip(axes, specs):
            ax.loglog(ells[2:],  cl_meas[2:],  '-',  color=color, alpha=0.6,
                      label='Measured (inv-var wtd)')
            ax.loglog(leff,      cl_h_meas,    'o',  color=color, alpha=0.9,
                      markersize=3, label='Measured (homogenized, full)')
            ax.loglog(leff,      cl_ha_meas,   's',  color='tab:purple', alpha=0.9,
                      markersize=3, label='Measured (homogenized, analysis mask)')
            ax.loglog(ells[2:],  Nl_th[2:],    '--', color='k',   alpha=0.8,
                      label='Theory (fit white + 1/f)')
            ax.set_ylabel(rf'$N_\ell^{{{label}}}$ [$\mu K^2$]')
            ax.set_title(f'{ch_name}  f={f:.2f} — {label}')
            ax.set_xlabel(r'$\ell$')
            ax.grid(True, which='both', alpha=0.2)
            ax.set_xlim([10, LMAX])
            ax.legend(fontsize='small')

        p_str = (f"lk_T={wp['ell_knee_i']}, α_T={wp['alpha_knee_i']}  |  "
                 f"lk_P={wp['ell_knee_p']}, α_P={wp['alpha_knee_p']}")
        plt.suptitle(
            f"{ch_name}  mixed wide+delens  f={f:.2f}  "
            f"(wide: {1-f:.2f}, delens: {f:.2f})\n"
            f"{p_str}\n"
            f"Meas depth: T={uKarcmin_TT:.1f}, P={uKarcmin_P:.1f} uK-arcmin   |   "
            f"Var-map expected: T={depth_exp_T:.1f}, P={depth_exp_P:.1f} uK-arcmin",
            fontsize=9, y=1.05,
        )
        plt.tight_layout()
        plt.savefig(
            os.path.join(validation_folder, f'power_spectrum_{ch_name}_f{f:.2f}.png'),
            dpi=150, bbox_inches='tight',
        )
        plt.close()

    print(f"\nDone. Validation plots saved to {validation_folder}")


if __name__ == '__main__':
    main()
