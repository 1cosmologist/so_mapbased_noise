#!/usr/bin/env python
"""
Validate noise simulations for LAT delensing survey by comparing measured power
spectra with theory.
"""

import os
import yaml
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import cmocean
import skytools as st
import pymaster as nmt
from tqdm import tqdm

# Configuration
NSIMS_FOR_AVG = 5
YAML_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/instr_params_lat-delens_channels.yaml'
BASE_OUTPUT_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/output'
RESOURCE_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources'
ANALYSIS_MASK_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/so_lat_deep_wide_analysis_2048.fits'
LMAX = 4500

# Expected white noise depths for reference lines (polarization sensitivity; T = P / sqrt(2))
EXPECTED_DEPTH_uKarcmin = {
    'LF027': 17.48577837911208,
    'LF039':  8.992686023543357,
    'MF093':  1.3006870704414162,
    'MF145':  1.400739922013833,
    'HF225':  3.325158341723625,
    'HF280':  8.34826987922102,
}


def get_output_folder(yaml_path, base_output_dir):
    """Extract output folder name from last two underscore-separated words of yaml filename."""
    yaml_basename = os.path.basename(yaml_path)
    yaml_name = os.path.splitext(yaml_basename)[0]
    parts = yaml_name.split('_')
    if len(parts) >= 2:
        folder_name = '_'.join(parts[-2:])
    else:
        folder_name = yaml_name
    return os.path.join(base_output_dir, folder_name)


def get_delens_files(ch_name, resource_dir):
    """Return relhits and binary mask paths for the delensing survey based on channel name."""
    if 'LF' in ch_name:
        band = 'lf'
    elif 'MF' in ch_name:
        band = 'mf'
    elif 'HF' in ch_name:
        band = 'uhf'
    else:
        raise ValueError(f"Unknown band in {ch_name}")

    relhits = os.path.join(resource_dir, f'so_lat_delens_wide_{band}_relhits_2048.fits')
    binary  = os.path.join(resource_dir, f'so_lat_delens_wide_{band}_bin_2048.fits')
    return relhits, binary


def uKarcmin2Nl(uKarcmin):
    """Convert noise level in uK-arcmin to Nl in uK^2-sr."""
    return (uKarcmin * np.deg2rad(1./60.))**2


def Nl2uKarcmin(Nl):
    """Convert Nl in uK^2-sr to uK-arcmin."""
    return np.sqrt(Nl) / np.deg2rad(1./60.)


def get_theory_Nl_ip(Nl_white_I, Nl_white_P, params, lmax):
    """
    Compute theory noise power spectrum with 1/f component for I and P.

    params: dict containing ell_knee_i, alpha_knee_i, ell_knee_p, alpha_knee_p
    """
    ells = np.arange(lmax + 1)

    factor_i = np.ones_like(ells, dtype=float)
    factor_i[2:] = 1. + (ells[2:] / params['ell_knee_i'])**params['alpha_knee_i']
    Nl_I = Nl_white_I * factor_i

    factor_p = np.ones_like(ells, dtype=float)
    factor_p[2:] = 1. + (ells[2:] / params['ell_knee_p'])**params['alpha_knee_p']
    Nl_P = Nl_white_P * factor_p

    return ells, Nl_I, Nl_P


def main():
    # Load channel parameters from YAML
    with open(YAML_FILE, 'r') as f:
        params = yaml.safe_load(f)

    channels = list(params.keys())
    print(f"Found {len(channels)} channels: {channels}")

    # Load analysis mask (applied to all noise maps and hit maps)
    print(f"Loading analysis mask: {ANALYSIS_MASK_FILE}")
    analysis_mask = hp.read_map(ANALYSIS_MASK_FILE)

    # Get input/output folders
    noise_folder = get_output_folder(YAML_FILE, BASE_OUTPUT_DIR)
    validation_folder = os.path.join(noise_folder, 'validation')
    os.makedirs(validation_folder, exist_ok=True)
    print(f"Noise maps folder:       {noise_folder}")
    print(f"Validation output folder: {validation_folder}")

    # Process each channel
    for ch_name in channels:
        print(f"\n{'='*60}")
        print(f"Processing channel: {ch_name}")
        print(f"{'='*60}")

        ch_params = params[ch_name]
        nside = ch_params.get('nside', 2048)

        # Get masking files (delensing-survey variants)
        relhits_file, binary_file = get_delens_files(ch_name, RESOURCE_DIR)
        print(f"  Using hits file: {os.path.basename(relhits_file)}")

        if not os.path.exists(relhits_file):
            print(f"  ERROR: Relhits file not found: {relhits_file}")
            continue

        # Load relative hits map for weighting; restrict to analysis mask
        relhits = hp.read_map(relhits_file) * analysis_mask
        weight = np.zeros_like(relhits)
        weight[relhits >= 1e-2] = relhits[relhits >= 1e-2]

        fsky = st.fsky(weight)

        binary_mask = hp.read_map(binary_file) * analysis_mask
        mask_apo = nmt.mask_apodization(binary_mask, 1.0, apotype='C2')

        print(f"  Effective f_sky after weighting: {fsky:.4f}")

        # =====================================================================
        # Part 1: Plot mollview of zeroth realization
        # =====================================================================
        print(f"  Creating mollview plots for realization 0...")

        noise_file_test = os.path.join(
            noise_folder,
            f"sobs_lat-noise_{ch_name}_mc000_nside{nside:04d}.fits"
        )
        if not os.path.exists(noise_file_test):
            print(f"  Noise file not found: {noise_file_test}")
            continue

        noise_map = hp.read_map(noise_file_test, field=None) * analysis_mask  # I, Q, U

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        titles = ['I (Temperature)', 'Q', 'U']

        for i, (ax, title) in enumerate(zip(axes, titles)):
            plt.sca(ax)
            map_data = noise_map[i]
            valid_pixels = map_data[binary_mask > 0.5]
            vmax = np.std(valid_pixels) * 3 if len(valid_pixels) > 0 else 100
            hp.mollview(
                map_data,
                title=f'{ch_name} - {title}',
                cmap=cmocean.cm.balance,
                min=-vmax, max=vmax,
                unit=r'$\mu K$',
                hold=True
            )
            hp.graticule(lw=0.5, color='gray', ls='--', alpha=0.5)

        plt.savefig(
            os.path.join(validation_folder, f'mollview_{ch_name}_mc000.png'),
            dpi=150, bbox_inches='tight'
        )
        plt.close()

        # =====================================================================
        # Part 2: Compute average power spectrum
        # =====================================================================
        print(f"  Computing average power spectrum over {NSIMS_FOR_AVG} realizations...")

        bins = nmt.NmtBin.from_lmax_linear(LMAX, 30)
        leff = bins.get_effective_ells()
        nmt_wsp = nmt.NmtWorkspace()

        cl_sum_TT = np.zeros(LMAX + 1)
        cl_sum_EE = np.zeros(LMAX + 1)
        cl_sum_BB = np.zeros(LMAX + 1)

        cl_h_sum_TT = np.zeros(len(leff))
        cl_h_sum_EE = np.zeros(len(leff))
        cl_h_sum_BB = np.zeros(len(leff))

        for sim_idx in tqdm(range(NSIMS_FOR_AVG), desc=f"  {ch_name} spectra", ncols=120):
            noise_file = os.path.join(
                noise_folder,
                f"sobs_lat-noise_{ch_name}_mc{sim_idx:03d}_nside{nside:04d}.fits"
            )
            if not os.path.exists(noise_file):
                print(f"  Missing sim {sim_idx}, skipping...")
                continue

            noise_map = hp.read_map(noise_file, field=None) * analysis_mask

            # Hits-weighted power spectra
            weighted_map = noise_map * weight
            cls = hp.anafast(weighted_map, lmax=LMAX, pol=True, nspec=3) / fsky

            cl_sum_TT += cls[0]
            cl_sum_EE += cls[1]
            cl_sum_BB += cls[2]

            # Homogenized noise for NaMaster decoupling
            homo_noise_map = noise_map * np.sqrt(weight)

            if sim_idx == 0:
                field = nmt.NmtField(
                    mask_apo, [homo_noise_map[0]],
                    purify_e=False, purify_b=False,
                    lmax=LMAX, lmax_mask=LMAX
                )
                nmt_wsp.compute_coupling_matrix(field, field, bins)

            cls_h = hp.anafast(homo_noise_map * mask_apo, lmax=LMAX, pol=True, nspec=3)

            cls_h_binned = np.zeros((3, len(leff)))
            for i in range(3):
                cls_h_binned[i] = nmt_wsp.decouple_cell(cls_h[i].reshape(1, -1))[0]

            cl_h_sum_TT += cls_h_binned[0]
            cl_h_sum_EE += cls_h_binned[1]
            cl_h_sum_BB += cls_h_binned[2]

        # Averages
        cl_avg_TT = cl_sum_TT / NSIMS_FOR_AVG
        cl_avg_EE = cl_sum_EE / NSIMS_FOR_AVG
        cl_avg_BB = cl_sum_BB / NSIMS_FOR_AVG

        cl_h_avg_TT = cl_h_sum_TT / NSIMS_FOR_AVG
        cl_h_avg_EE = cl_h_sum_EE / NSIMS_FOR_AVG
        cl_h_avg_BB = cl_h_sum_BB / NSIMS_FOR_AVG

        # =====================================================================
        # Part 3: Estimate white noise level and compute theory curves
        # =====================================================================
        print(f"  Estimating white noise level and computing theory...")

        idx_min = 3000
        idx_max = min(4500, LMAX)

        Nl_white_TT_meas = np.mean(cl_avg_TT[idx_min:idx_max])
        Nl_white_P_meas  = np.mean([
            np.mean(cl_avg_EE[idx_min:idx_max]),
            np.mean(cl_avg_BB[idx_min:idx_max])
        ])

        uKarcmin_TT = Nl2uKarcmin(Nl_white_TT_meas)
        uKarcmin_P  = Nl2uKarcmin(Nl_white_P_meas)

        print(f"  Estimated Depth (from spectra ell={idx_min}-{idx_max}):")
        print(f"    T: {uKarcmin_TT:.2f} uK-arcmin")
        print(f"    P: {uKarcmin_P:.2f} uK-arcmin")
        print(f"    Ratio P/T: {uKarcmin_P/uKarcmin_TT:.2f} (Exp: ~1.414)")

        ells, Nl_theory_T, Nl_theory_P = get_theory_Nl_ip(
            Nl_white_TT_meas, Nl_white_P_meas, ch_params, LMAX
        )

        # =====================================================================
        # Part 4: Plot
        # =====================================================================
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        p_str = (f"lk_T={ch_params['ell_knee_i']}, a_T={ch_params['alpha_knee_i']} | "
                 f"lk_P={ch_params['ell_knee_p']}, a_P={ch_params['alpha_knee_p']}")

        # Reference horizontal lines from expected depths
        ch_prefix = ch_name[:5]
        Nl_ref_T = Nl_ref_P = None
        if ch_prefix in EXPECTED_DEPTH_uKarcmin:
            depth_P = EXPECTED_DEPTH_uKarcmin[ch_prefix]
            depth_T = depth_P / np.sqrt(2)
            Nl_ref_T = uKarcmin2Nl(depth_T)
            Nl_ref_P = uKarcmin2Nl(depth_P)

        # TT
        ax = axes[0]
        ax.loglog(ells[2:], cl_avg_TT[2:], '-',  label='Measured (hits wtd)',    alpha=0.6, color='tab:blue')
        ax.loglog(leff,     cl_h_avg_TT,   'o',  label='Measured (homogenized)', alpha=0.9, color='tab:blue', markersize=3)
        ax.loglog(ells[2:], Nl_theory_T[2:], '--', label='Theory (Fit White)',   alpha=0.8, color='k')
        if Nl_ref_T is not None:
            ax.axhline(Nl_ref_T, ls=':', color='red', alpha=0.8,
                       label=f'Expected ({depth_T:.2f} uK-am)')
        ax.set_ylabel(r'$N_\ell^{TT}$ [$\mu K^2$]')
        ax.set_title(f'{ch_name} - TT')
        ax.legend(fontsize='small')

        # EE
        ax = axes[1]
        ax.loglog(ells[2:], cl_avg_EE[2:], '-',  label='Meas',     alpha=0.6, color='tab:orange')
        ax.loglog(leff,     cl_h_avg_EE,   'o',  label='Meas (hom)', alpha=0.9, color='tab:orange', markersize=3)
        ax.loglog(ells[2:], Nl_theory_P[2:], '--', label='Theory',  alpha=0.8, color='k')
        if Nl_ref_P is not None:
            ax.axhline(Nl_ref_P, ls=':', color='red', alpha=0.8,
                       label=f'Expected ({depth_P:.2f} uK-am)')
        ax.set_ylabel(r'$N_\ell^{EE}$ [$\mu K^2$]')
        ax.set_title(f'{ch_name} - EE')

        # BB
        ax = axes[2]
        ax.loglog(ells[2:], cl_avg_BB[2:], '-',  label='Meas',     alpha=0.6, color='tab:green')
        ax.loglog(leff,     cl_h_avg_BB,   'o',  label='Meas (hom)', alpha=0.9, color='tab:green', markersize=3)
        ax.loglog(ells[2:], Nl_theory_P[2:], '--', label='Theory',  alpha=0.8, color='k')
        if Nl_ref_P is not None:
            ax.axhline(Nl_ref_P, ls=':', color='red', alpha=0.8,
                       label=f'Expected ({depth_P:.2f} uK-am)')
        ax.set_ylabel(r'$N_\ell^{BB}$ [$\mu K^2$]')
        ax.set_title(f'{ch_name} - BB')

        for ax in axes:
            ax.set_xlabel(r'$\ell$')
            ax.grid(True, which='both', alpha=0.2)
            ax.set_xlim([10, LMAX])
            ax.legend(fontsize='small')

        plt.suptitle(
            f'{ch_name} Delensing Noise Validation\n{p_str}\n'
            f'Est Depths: T={uKarcmin_TT:.1f}, P={uKarcmin_P:.1f} uK-arcmin',
            fontsize=10, y=1.05
        )
        plt.tight_layout()
        plt.savefig(
            os.path.join(validation_folder, f'power_spectrum_{ch_name}.png'),
            dpi=150, bbox_inches='tight'
        )
        plt.close()

    print(f"\nValidation complete! Results saved to: {validation_folder}")


if __name__ == '__main__':
    main()
