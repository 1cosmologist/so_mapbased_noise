#!/usr/bin/env python
"""
Validate noise simulations by comparing measured power spectra with theory.
"""

import os
import yaml
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import cmocean
import skytools as st
from tqdm import tqdm

# Configuration
NSIMS_FOR_AVG = 5
YAML_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/instr_params_baseline_pessimistic.yaml'
BASE_OUTPUT_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/output'
RELHITS_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/so_sat_relhits_C_nside512.fits'
FULL_BINARY_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/so_sat_full-binary_C_nside512.fits'
LMAX = 1000


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


def uKarcmin2Nl(uKarcmin):
    """Convert noise level in uK-arcmin to Nl in uK^2-sr."""
    return (uKarcmin * np.deg2rad(1./60.))**2


def get_theory_Nl(uKarcmin, ell_knee, alpha_knee, lmax):
    """Compute theory noise power spectrum with 1/f component."""
    ells = np.arange(lmax + 1)
    Nl = np.zeros(lmax + 1)
    Nl[2:] = uKarcmin2Nl(uKarcmin) * (1. + (ells[2:] / ell_knee)**alpha_knee)
    return ells, Nl


def main():
    # Load channel parameters from YAML
    with open(YAML_FILE, 'r') as f:
        params = yaml.safe_load(f)
    
    channels = list(params.keys())
    print(f"Found {len(channels)} channels: {channels}")
    
    # Get input/output folders
    noise_folder = get_output_folder(YAML_FILE, BASE_OUTPUT_DIR)
    validation_folder = os.path.join(noise_folder, 'validation')
    os.makedirs(validation_folder, exist_ok=True)
    print(f"Noise maps folder: {noise_folder}")
    print(f"Validation output folder: {validation_folder}")
    
    # Load relative hits map for weighting
    relhits = hp.read_map(RELHITS_FILE)
    # Create weight map (zero where hits are too low)
    weight = np.zeros_like(relhits)
    weight[relhits >= 1e-2] = relhits[relhits >= 1e-2]
    fsky = st.fsky(weight)
    
    binary_mask = hp.read_map(FULL_BINARY_FILE)
    fsky_binary = st.fsky(binary_mask)
    
    print(f"Effective f_sky after weighting: {fsky:.4f}")
    print(f"Effective f_sky for binary mask: {fsky_binary:.4f}")
    
    # Process each channel
    for ch_name in channels:
        print(f"\n{'='*60}")
        print(f"Processing channel: {ch_name}")
        print(f"{'='*60}")
        
        ch_params = params[ch_name]
        nside = ch_params['nside']
        uKarcmin = ch_params['noise_uK_arcmin']
        ell_knee = ch_params['ell_knee']
        alpha_knee = ch_params['alpha_knee']
        
        # =====================================================================
        # Part 1: Plot mollview of zeroth realization
        # =====================================================================
        print(f"  Creating mollview plots for realization 0...")
        
        noise_file = os.path.join(
            noise_folder,
            f"sobs_noise_{ch_name}_mc000_nside{nside:04d}.fits"
        )
        noise_map = hp.read_map(noise_file, field=None)  # Read all 3 fields (I, Q, U)
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        titles = ['I (Temperature)', 'Q', 'U']
        
        for i, (ax, title) in enumerate(zip(axes, titles)):
            plt.sca(ax)
            # Get map range for symmetric colorbar
            map_data = noise_map[i]
            vmax = np.std(map_data[map_data != 0]) * 3
            
            hp.mollview(
                map_data,
                title=f'{ch_name} - {title} (mc000)',
                cmap=cmocean.cm.balance,
                min=-vmax, max=vmax,
                unit=r'$\mu K$',
                hold=True
            )
            
            hp.graticule(lw=0.5, color='gray', ls='--', alpha=0.5)
        
        # plt.tight_layout()
        plt.savefig(os.path.join(validation_folder, f'mollview_{ch_name}_mc000.png'), 
                    dpi=150, bbox_inches='tight')
        plt.close()
        
        # =====================================================================
        # Part 2: Compute average power spectrum over first 25 sims
        # =====================================================================
        print(f"  Computing average power spectrum over {NSIMS_FOR_AVG} realizations...")
        
        # Initialize arrays for accumulating spectra (TT, EE, BB, TE, EB, TB)
        cl_sum_TT = np.zeros(LMAX + 1)
        cl_sum_EE = np.zeros(LMAX + 1)
        cl_sum_BB = np.zeros(LMAX + 1)
        
        cl_h_sum_TT = np.zeros(LMAX + 1)
        cl_h_sum_EE = np.zeros(LMAX + 1)
        cl_h_sum_BB = np.zeros(LMAX + 1)
        
        for sim_idx in tqdm(range(NSIMS_FOR_AVG), desc=f"  {ch_name} spectra", ncols=120):
            noise_file = os.path.join(
                noise_folder,
                f"sobs_noise_{ch_name}_mc{sim_idx:03d}_nside{nside:04d}.fits"
            )
            noise_map = hp.read_map(noise_file, field=None)
            
            # Apply relhits weighting
            weighted_map = noise_map * weight
            
            # Compute power spectrum with anafast
            cls = hp.anafast(weighted_map, lmax=LMAX, pol=True, nspec=3) / fsky
            
            homo_noise_map = noise_map * np.sqrt(weight)
            cls_h = hp.anafast(homo_noise_map, lmax=LMAX, pol=True, nspec=3) / fsky_binary
            # cls returns: TT, EE, BB
            
            cl_sum_TT += cls[0]
            cl_sum_EE += cls[1]
            cl_sum_BB += cls[2]
            
            cl_h_sum_TT += cls_h[0]
            cl_h_sum_EE += cls_h[1]
            cl_h_sum_BB += cls_h[2]
        
        # Average spectra
        cl_avg_TT = cl_sum_TT / NSIMS_FOR_AVG
        cl_avg_EE = cl_sum_EE / NSIMS_FOR_AVG
        cl_avg_BB = cl_sum_BB / NSIMS_FOR_AVG
        
        cl_h_avg_TT = cl_h_sum_TT / NSIMS_FOR_AVG
        cl_h_avg_EE = cl_h_sum_EE / NSIMS_FOR_AVG
        cl_h_avg_BB = cl_h_sum_BB / NSIMS_FOR_AVG
        
        # =====================================================================
        # Part 3: Compute theory spectrum and bin both
        # =====================================================================
        print(f"  Computing theory spectrum and binning...")
        
        ells = np.arange(LMAX + 1)
        ells_theory, Nl_theory = get_theory_Nl(uKarcmin, ell_knee, alpha_knee, LMAX) # / np.sqrt(0.85)
        
        # Nl for polarization is 2x temperature (factor of 2 for Q and U)
        Nl_theory_T = Nl_theory / 2  # Temperature gets factor of 1/sqrt(2) in code
        Nl_theory_P = Nl_theory      # Polarization
        
        # =====================================================================
        # Part 4: Plot comparison
        # =====================================================================
        print(f"  Creating power spectrum comparison plots...")
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        
        # TT spectrum
        ax = axes[0]
        ax.semilogy(ells[2:], cl_avg_TT[2:], '-', label='Measured (hits weighted)', alpha=0.8)
        ax.semilogy(ells[2:], cl_h_avg_TT[2:], '--', label='Measured (homogenized)', alpha=0.8)
        ax.semilogy(ells_theory[2:], Nl_theory_T[2:], '-', label='Theory', alpha=0.7)
        ax.set_xlabel(r'$\ell$')
        ax.set_ylabel(r'$C_\ell^{TT}$ [$\mu K^2$]')
        ax.set_title(f'{ch_name} - TT Power Spectrum')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0, LMAX])
        
        # EE spectrum
        ax = axes[1]
        ax.semilogy(ells[2:], cl_avg_EE[2:], '-', label='Measured (hits weighted)', alpha=0.8)
        ax.semilogy(ells[2:], cl_h_avg_EE[2:], '--', label='Measured (homogenized)', alpha=0.8)
        ax.semilogy(ells_theory[2:], Nl_theory_P[2:], '-', label='Theory', alpha=0.7)
        ax.set_xlabel(r'$\ell$')
        ax.set_ylabel(r'$C_\ell^{EE}$ [$\mu K^2$]')
        ax.set_title(f'{ch_name} - EE Power Spectrum')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0, LMAX])
        
        # BB spectrum
        ax = axes[2]
        ax.semilogy(ells[2:], cl_avg_BB[2:], '-', label='Measured (hits weighted)', alpha=0.8)
        ax.semilogy(ells[2:], cl_h_avg_BB[2:], '--', label='Measured (homogenized)', alpha=0.8)
        ax.semilogy(ells_theory[2:], Nl_theory_P[2:], '-', label='Theory', alpha=0.7)
        ax.set_xlabel(r'$\ell$')
        ax.set_ylabel(r'$C_\ell^{BB}$ [$\mu K^2$]')
        ax.set_title(f'{ch_name} - BB Power Spectrum')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0, LMAX])
        
        plt.suptitle(f'{ch_name}: Noise = {uKarcmin} uK-arcmin, ell_knee = {ell_knee}, alpha = {alpha_knee}',
                     fontsize=12, y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(validation_folder, f'power_spectrum_{ch_name}.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"  Done with {ch_name}")
    
    print(f"\n{'='*60}")
    print(f"Validation complete! Results saved to: {validation_folder}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
