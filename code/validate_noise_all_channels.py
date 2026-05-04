#!/usr/bin/env python
"""
Validate noise simulations for all channels (variance_map method).
Compares measured power spectra with fitted noise model.
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
from scipy.optimize import curve_fit

# Configuration
NSIMS = 50
YAML_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/instr_params_all_channels.yaml'
NOISE_MAP_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/output/all_channels'
FIGURES_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/figures'
RELHITS_SOUTH_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/so_south_relhits_2dC2.fits'
SOUTH_BINARY_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/so_south_patch_binary.fits'
LMAX = 1000


def noise_model(ells, N0, ell_knee, alpha):
    """
    Noise power spectrum model: N_ell = N0 * [1 + (ell/ell_knee)^alpha]
    
    Parameters
    ----------
    ells : array
        Multipole values
    N0 : float
        White noise level (in power spectrum units)
    ell_knee : float
        Knee multipole where atmospheric noise equals white noise
    alpha : float
        Spectral slope of atmospheric component (should be negative)
    
    Returns
    -------
    array
        Model noise power spectrum
    """
    return N0 * (1.0 + (ells / ell_knee)**alpha)


def fit_N0_only(ells, cl_measured):
    """
    Fit white noise level N0 only (assumes high ell where 1/f is negligible).
    
    At high ell, the model reduces to N_ell ~ N0, so we just fit a constant.
    
    Parameters
    ----------
    ells : array
        Multipole values (should be high ell where noise is white)
    cl_measured : array
        Measured power spectrum
    
    Returns
    -------
    float
        Best-fit N0 value
    """
    # Simple mean of the spectrum at high ell where it should be flat
    N0_fit = np.mean(cl_measured)
    return N0_fit


def fit_ell_knee_alpha(ells, cl_measured, N0_fixed, p0=None):
    """
    Fit ell_knee and alpha with N0 fixed.
    
    Model: N_ell = N0 * [1 + (ell/ell_knee)^alpha]
    
    Parameters
    ----------
    ells : array
        Multipole values
    cl_measured : array
        Measured power spectrum
    N0_fixed : float
        Fixed white noise level from high-ell fit
    p0 : tuple, optional
        Initial guess for (ell_knee, alpha)
    
    Returns
    -------
    tuple
        Best-fit parameters (ell_knee, alpha) and covariance matrix
    """
    def model_fixed_N0(ells, ell_knee, alpha):
        return N0_fixed * (1.0 + (ells / ell_knee)**alpha)
    
    if p0 is None:
        p0 = (50.0, -3.0)
    
    # Bounds: ell_knee > 1, alpha < 0
    bounds = ([1, -10], [500, 0])
    
    try:
        popt, pcov = curve_fit(model_fixed_N0, ells, cl_measured, p0=p0, bounds=bounds, maxfev=5000)
        return popt, pcov
    except Exception as e:
        print(f"    Warning: Fit failed with error: {e}")
        return None, None


def uKarcmin2Nl(uKarcmin):
    """Convert noise level in uK-arcmin to Nl in uK^2-sr."""
    return (uKarcmin * np.deg2rad(1./60.))**2


def main():
    # Create figures directory
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(os.path.join(FIGURES_DIR, 'all_channels_validation'), exist_ok=True)
    validation_dir = os.path.join(FIGURES_DIR, 'all_channels_validation')
    
    # Load channel parameters from YAML
    with open(YAML_FILE, 'r') as f:
        params = yaml.safe_load(f)
    
    channels = list(params.keys())
    print(f"Found {len(channels)} channels: {channels}")
    print(f"Noise maps directory: {NOISE_MAP_DIR}")
    print(f"Figures output directory: {validation_dir}")
    
    # Load south patch relative hits and binary mask
    relhits_south = hp.read_map(RELHITS_SOUTH_FILE)
    binary_south = hp.read_map(SOUTH_BINARY_FILE)
    
    # Create weight map for south patch
    weight_south = relhits_south.copy()
    fsky_south = st.fsky(weight_south)
    
    # Apodized mask for NaMaster
    mask_apo = nmt.mask_apodization(binary_south, 2., apotype='C2')
    
    print(f"South patch f_sky (weighted): {fsky_south:.4f}")
    
    # Store fit results for all channels
    fit_results = {}
    
    # Process each channel
    for ch_name in channels:
        print(f"\n{'='*60}")
        print(f"Processing channel: {ch_name}")
        print(f"{'='*60}")
        
        ch_params = params[ch_name]
        nside = ch_params['nside']
        ell_knee_config = ch_params['ell_knee']
        alpha_config = ch_params['alpha_knee']
        tel_yrs = ch_params['tel_yrs']
        
        # Build filename pattern for this channel
        # Format: sobs_noise_{ch_name}_mc{sim_idx:03d}_nside{nside:04d}.fits
        
        # =====================================================================
        # Part 1: Plot mollview of zeroth realization (I, Q, U as subplots)
        # =====================================================================
        print(f"  Creating mollview plots for realization 0...")
        
        noise_file = os.path.join(
            NOISE_MAP_DIR,
            f"sobs_noise_{ch_name}_mc000_nside{nside:04d}.fits"
        )
        
        if not os.path.exists(noise_file):
            print(f"  Warning: File not found: {noise_file}")
            print(f"  Skipping channel {ch_name}")
            continue
        
        noise_map = hp.read_map(noise_file, field=None)  # Read all 3 fields (I, Q, U)
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        titles = ['I (Temperature)', 'Q', 'U']
        
        for i, (ax, title) in enumerate(zip(axes, titles)):
            plt.sca(ax)
            # Get 2-sigma range for symmetric colorbar (use pixels where relhits > 0)
            map_data = noise_map[i]
            valid_data = map_data[relhits_south > 0]
            sigma = np.std(valid_data)
            vmax = 2 * sigma
            
            hp.mollview(
                map_data,
                title=f'{ch_name} - {title} (mc000)',
                cmap=cmocean.cm.balance,
                min=-vmax, max=vmax,
                unit=r'$\mu K$',
                hold=True
            )
            hp.graticule(lw=0.5, color='gray', ls='--', alpha=0.5)
        
        plt.suptitle(f'{ch_name}: {tel_yrs:.2f} tel-yrs, ell_knee={ell_knee_config}, alpha={alpha_config}', 
                     fontsize=12, y=1.02)
        plt.savefig(os.path.join(validation_dir, f'mollview_{ch_name}_mc000.png'), 
                    dpi=150, bbox_inches='tight')
        plt.close()
        
        # =====================================================================
        # Part 2: Compute average BB power spectrum over all sims (south patch)
        # =====================================================================
        print(f"  Computing average BB power spectrum over {NSIMS} realizations...")
        
        # Set up NaMaster binning
        bins = nmt.NmtBin.from_lmax_linear(LMAX, 15)
        leff = bins.get_effective_ells()
        nmt_wsp_relhits = nmt.NmtWorkspace()  # workspace for relhits-weighted (uses weight_south as mask)
        nmt_wsp_homo = nmt.NmtWorkspace()      # workspace for homogenized (uses mask_apo)
        
        # Arrays for accumulating spectra (4 methods: anafast/NaMaster x relhits/homogenized)
        cl_sum_BB_anafast_relhits = np.zeros(LMAX + 1)   # relhits-weighted anafast
        cl_sum_BB_anafast_homo = np.zeros(LMAX + 1)      # homogenized anafast
        cl_sum_BB_nmt_relhits = np.zeros(len(leff))       # relhits-weighted NaMaster (decoupled)
        cl_sum_BB_nmt_homo = np.zeros(len(leff))          # homogenized NaMaster (decoupled)
        
        ells = np.arange(LMAX + 1)
        fsky_homo = st.fsky(mask_apo)
        
        for sim_idx in tqdm(range(NSIMS), desc=f"  {ch_name} spectra", ncols=120):
            noise_file = os.path.join(
                NOISE_MAP_DIR,
                f"sobs_noise_{ch_name}_mc{sim_idx:03d}_nside{nside:04d}.fits"
            )
            noise_map = hp.read_map(noise_file, field=None)
            
            # Relhits-weighted map (uses weight_south as both weighting and effective mask)
            weighted_map = noise_map * weight_south
            
            # Homogenized map (uses sqrt(weight) to homogenize, then mask_apo for geometry)
            homo_noise_map = noise_map * np.sqrt(weight_south)
            
            # Method 1: relhits-weighted anafast
            cls_anafast_relhits = hp.anafast(weighted_map, lmax=LMAX, pol=True)
            # cls_anafast returns: TT, EE, BB, TE, EB, TB for pol=True
            cl_sum_BB_anafast_relhits += cls_anafast_relhits[2] / fsky_south
            
            # Method 2: homogenized anafast (apply mask_apo before anafast)
            cls_anafast_homo = hp.anafast(homo_noise_map * mask_apo, lmax=LMAX, pol=True)
            cl_sum_BB_anafast_homo += cls_anafast_homo[2] / fsky_homo
            
            if sim_idx == 0:
                # Compute coupling matrix for relhits-weighted case (uses weight_south as mask)
                field_relhits = nmt.NmtField(weight_south, [noise_map[0]], 
                                             purify_e=False, purify_b=False, lmax=LMAX, lmax_mask=LMAX)
                nmt_wsp_relhits.compute_coupling_matrix(field_relhits, field_relhits, bins)
                
                # Compute coupling matrix for homogenized case (uses mask_apo)
                field_homo = nmt.NmtField(mask_apo, [homo_noise_map[0]], 
                                          purify_e=False, purify_b=False, lmax=LMAX, lmax_mask=LMAX)
                nmt_wsp_homo.compute_coupling_matrix(field_homo, field_homo, bins)
            
            # Method 3: relhits-weighted NaMaster
            # Use anafast to get pseudo-Cl (weighted_map already includes weight_south), then decouple
            cls_relhits_pseudo = hp.anafast(weighted_map, lmax=LMAX, pol=True)
            cls_relhits_binned = nmt_wsp_relhits.decouple_cell(cls_relhits_pseudo[2].reshape(1, -1))[0]
            cl_sum_BB_nmt_relhits += cls_relhits_binned
            
            # Method 4: homogenized NaMaster
            # Use anafast to get pseudo-Cl (homo_noise_map * mask_apo), then decouple
            cls_homo_pseudo = hp.anafast(homo_noise_map * mask_apo, lmax=LMAX, pol=True)
            cls_homo_binned = nmt_wsp_homo.decouple_cell(cls_homo_pseudo[2].reshape(1, -1))[0]
            cl_sum_BB_nmt_homo += cls_homo_binned
        
        # Average spectra
        cl_avg_BB_anafast_relhits = cl_sum_BB_anafast_relhits / NSIMS
        cl_avg_BB_anafast_homo = cl_sum_BB_anafast_homo / NSIMS
        cl_avg_BB_nmt_relhits = cl_sum_BB_nmt_relhits / NSIMS
        cl_avg_BB_nmt_homo = cl_sum_BB_nmt_homo / NSIMS
        
        # =====================================================================
        # Part 3: Two-step fit of noise model to homogenized spectrum
        # =====================================================================
        print(f"  Fitting noise model N_ell = N_0[1 + (ell/ell_knee)^alpha]...")
        
        # Step 1: Fit N0 from high ell (500-1000) where noise is white
        high_ell_mask = (leff > 500) & (leff < 1000)
        if np.sum(high_ell_mask) > 0:
            N0_fit = fit_N0_only(leff[high_ell_mask], cl_avg_BB_nmt_homo[high_ell_mask])
            print(f"  Step 1: Fit N_0 from ell=[500, 1000]")
            uKarcmin_fit = np.sqrt(N0_fit) / np.deg2rad(1./60.)
            print(f"    N_0 = {N0_fit:.4e} uK^2-sr (equiv. {uKarcmin_fit:.2f} uK-arcmin)")
        else:
            print(f"  Warning: Not enough bins in ell=[500, 1000] for N0 fit")
            N0_fit = None
        
        # Use ell_knee and alpha from config file (not fitted)
        if N0_fit is not None:
            print(f"  Using config values: ell_knee = {ell_knee_config}, alpha = {alpha_config}")
            print(f"  Fit results:")
            print(f"    N_0 = {N0_fit:.4e} uK^2-sr (equiv. {uKarcmin_fit:.2f} uK-arcmin)")
            print(f"    ell_knee = {ell_knee_config} (from config)")
            print(f"    alpha = {alpha_config} (from config)")
            
            fit_results[ch_name] = {
                'N0': N0_fit,
                'uKarcmin_equiv': uKarcmin_fit,
                'ell_knee': ell_knee_config,
                'alpha': alpha_config,
                'ell_knee_config': ell_knee_config,
                'alpha_config': alpha_config
            }
            
            # Generate model curve using config values for ell_knee and alpha
            ells_model = np.arange(2, LMAX + 1)
            model_curve = noise_model(ells_model, N0_fit, ell_knee_config, alpha_config)
        else:
            fit_results[ch_name] = None
            model_curve = None
        
        # =====================================================================
        # Part 4: Plot BB power spectrum comparison
        # =====================================================================
        print(f"  Creating BB power spectrum plot...")
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Relhits-weighted case (blue): anafast line + NaMaster markers
        ax.loglog(ells[2:], cl_avg_BB_anafast_relhits[2:], '-', alpha=0.6, 
                  label='Relhits-weighted (anafast)', color='C0')
        ax.loglog(leff, cl_avg_BB_nmt_relhits, 'o', markersize=3, alpha=0.8,
                  label='Relhits-weighted (NaMaster)', color='C0')
        
        # Homogenized case (orange): anafast line + NaMaster markers
        ax.loglog(ells[2:], cl_avg_BB_anafast_homo[2:], '-', alpha=0.6, 
                  label='Homogenized (anafast)', color='C1')
        ax.loglog(leff, cl_avg_BB_nmt_homo, 'o', markersize=3, alpha=0.8,
                  label='Homogenized (NaMaster)', color='C1')
        
        # Plot noise model (N_0 fitted, ell_knee and alpha from config)
        if model_curve is not None:
            ax.loglog(ells_model, model_curve, '--', lw=2, 
                      label=f'Model: N₀={N0_fit:.2e}, ℓ_knee={ell_knee_config} (cfg), α={alpha_config} (cfg)', 
                      color='C2')
        
        ax.set_xlabel(r'$\ell$', fontsize=12)
        ax.set_ylabel(r'$C_\ell^{BB}$ [$\mu K^2$]', fontsize=12)
        ax.set_title(f'{ch_name} - BB Noise Power Spectrum (South Patch)\n'
                     f'{tel_yrs:.2f} tel-yrs, Config: ℓ_knee={ell_knee_config}, α={alpha_config}',
                     fontsize=11)
        ax.legend(fontsize=9, ncol=2)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([10, LMAX])
        
        plt.tight_layout()
        plt.savefig(os.path.join(validation_dir, f'bb_spectrum_{ch_name}.png'),
                    dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"  Done with {ch_name}")
    
    # =========================================================================
    # Summary of fit results
    # =========================================================================
    print(f"\n{'='*60}")
    print("SUMMARY OF RESULTS")
    print("(N_0 fitted from ell=500-1000; ell_knee and alpha from config)")
    print(f"{'='*60}")
    print(f"{'Channel':<35} {'N0 (uK-arcmin)':<18} {'ell_knee (cfg)':<18} {'alpha (cfg)':<15}")
    print("-" * 86)
    
    for ch_name, result in fit_results.items():
        if result is not None:
            print(f"{ch_name:<35} {result['uKarcmin_equiv']:<18.2f} {result['ell_knee']:<18} {result['alpha']:<15}")
        else:
            print(f"{ch_name:<35} {'FAILED':<18} {'-':<18} {'-':<15}")
    
    # =========================================================================
    # N_0 ratio comparison for frequency pairs with same tel_yrs
    # The ratio of N_0 should equal the ratio of NET^2
    # =========================================================================
    print(f"\n{'='*60}")
    print("N_0 RATIO COMPARISON (same tel_yrs pairs)")
    print("Expected: N_0 ratio = NET^2 ratio")
    print(f"{'='*60}")
    
    # NET values from SO specs (in uK*sqrt(s))
    NET_telescope = {27: 21, 39: 13, 93: 3.4, 145: 5.0, 225: 8.6, 280: 22, 346: 50}
    
    # Group channels by tel_yrs
    tel_yrs_groups = {}
    for ch_name, result in fit_results.items():
        if result is not None:
            tel_yrs = params[ch_name]['tel_yrs']
            freq = params[ch_name]['central_freq_GHz']
            if tel_yrs not in tel_yrs_groups:
                tel_yrs_groups[tel_yrs] = []
            tel_yrs_groups[tel_yrs].append((ch_name, freq, result['N0']))
    
    # For each tel_yrs group with multiple frequencies, compute N_0 ratios
    print(f"\n{'Pair':<50} {'tel_yrs':<12} {'N0 ratio (meas)':<18} {'NET^2 ratio (exp)':<18} {'Match?':<10}")
    print("-" * 108)
    
    for tel_yrs, channels in sorted(tel_yrs_groups.items()):
        if len(channels) >= 2:
            # Compare all pairs within this tel_yrs group
            for i in range(len(channels)):
                for j in range(i+1, len(channels)):
                    ch1_name, freq1, N0_1 = channels[i]
                    ch2_name, freq2, N0_2 = channels[j]
                    
                    # Compute measured N_0 ratio
                    N0_ratio_meas = N0_1 / N0_2
                    
                    # Compute expected NET^2 ratio
                    NET1 = NET_telescope.get(int(freq1), None)
                    NET2 = NET_telescope.get(int(freq2), None)
                    
                    if NET1 is not None and NET2 is not None:
                        NET2_ratio_exp = (NET1 / NET2)**2
                        
                        # Check if they match (within 10%)
                        ratio_diff = abs(N0_ratio_meas - NET2_ratio_exp) / NET2_ratio_exp * 100
                        match = "YES" if ratio_diff < 10 else f"NO ({ratio_diff:.1f}%)"
                        
                        pair_str = f"{ch1_name} / {ch2_name}"
                        print(f"{pair_str:<50} {tel_yrs:<12.2f} {N0_ratio_meas:<18.4f} {NET2_ratio_exp:<18.4f} {match:<10}")
    
    # =========================================================================
    # N_0 ratio comparison for same frequency with different tel_yrs
    # The ratio of N_0 should equal the inverse ratio of tel_yrs
    # =========================================================================
    print(f"\n{'='*60}")
    print("N_0 RATIO COMPARISON (same frequency, different tel_yrs)")
    print("Expected: N_0 ratio = (tel_yrs_2 / tel_yrs_1) = 1 / (tel_yrs ratio)")
    print(f"{'='*60}")
    
    # Group channels by frequency
    freq_groups = {}
    for ch_name, result in fit_results.items():
        if result is not None:
            tel_yrs = params[ch_name]['tel_yrs']
            freq = params[ch_name]['central_freq_GHz']
            if freq not in freq_groups:
                freq_groups[freq] = []
            freq_groups[freq].append((ch_name, tel_yrs, result['N0']))
    
    # For each frequency group with multiple tel_yrs, compute N_0 ratios
    print(f"\n{'Pair':<50} {'Freq (GHz)':<12} {'N0 ratio (meas)':<18} {'1/tel_yrs ratio (exp)':<22} {'Match?':<10}")
    print("-" * 112)
    
    for freq, channels in sorted(freq_groups.items()):
        if len(channels) >= 2:
            # Sort by tel_yrs for consistent ordering
            channels_sorted = sorted(channels, key=lambda x: x[1])
            
            # Compare all pairs within this frequency group
            for i in range(len(channels_sorted)):
                for j in range(i+1, len(channels_sorted)):
                    ch1_name, tel_yrs_1, N0_1 = channels_sorted[i]
                    ch2_name, tel_yrs_2, N0_2 = channels_sorted[j]
                    
                    # Compute measured N_0 ratio
                    N0_ratio_meas = N0_1 / N0_2
                    
                    # Compute expected inverse tel_yrs ratio
                    # N_0 ~ 1/t_obs, so N0_1/N0_2 = t_obs_2/t_obs_1 = tel_yrs_2/tel_yrs_1
                    inv_telyrs_ratio_exp = tel_yrs_2 / tel_yrs_1
                    
                    # Check if they match (within 10%)
                    ratio_diff = abs(N0_ratio_meas - inv_telyrs_ratio_exp) / inv_telyrs_ratio_exp * 100
                    match = "YES" if ratio_diff < 10 else f"NO ({ratio_diff:.1f}%)"
                    
                    pair_str = f"{ch1_name} / {ch2_name}"
                    print(f"{pair_str:<50} {freq:<12.0f} {N0_ratio_meas:<18.4f} {inv_telyrs_ratio_exp:<22.4f} {match:<10}")
    
    # =========================================================================
    # Save results to text file
    # =========================================================================
    results_file = os.path.join(validation_dir, 'fit_results.txt')
    with open(results_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("SO SAT Noise Validation Results\n")
        f.write("=" * 80 + "\n\n")
        
        # Noise model parameters (N_0 fitted, ell_knee and alpha from config)
        f.write("NOISE MODEL PARAMETERS\n")
        f.write("(N_0 fitted from ell=500-1000; ell_knee and alpha from config)\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Channel':<35} {'N0 (uK-arcmin)':<18} {'ell_knee (cfg)':<18} {'alpha (cfg)':<15}\n")
        f.write("-" * 86 + "\n")
        
        for ch_name, result in fit_results.items():
            if result is not None:
                f.write(f"{ch_name:<35} {result['uKarcmin_equiv']:<18.2f} {result['ell_knee']:<18} {result['alpha']:<15}\n")
            else:
                f.write(f"{ch_name:<35} {'FAILED':<18} {'-':<18} {'-':<15}\n")
        
        # N_0 ratio test: same tel_yrs, different frequencies
        f.write("\n\n" + "=" * 80 + "\n")
        f.write("N_0 RATIO TEST: Same tel_yrs, different frequencies\n")
        f.write("Expected: N_0 ratio = NET^2 ratio\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Pair':<50} {'tel_yrs':<12} {'N0 ratio (meas)':<18} {'NET^2 ratio (exp)':<18} {'Match?':<10}\n")
        f.write("-" * 108 + "\n")
        
        for tel_yrs, channels in sorted(tel_yrs_groups.items()):
            if len(channels) >= 2:
                for i in range(len(channels)):
                    for j in range(i+1, len(channels)):
                        ch1_name, freq1, N0_1 = channels[i]
                        ch2_name, freq2, N0_2 = channels[j]
                        
                        N0_ratio_meas = N0_1 / N0_2
                        NET1 = NET_telescope.get(int(freq1), None)
                        NET2 = NET_telescope.get(int(freq2), None)
                        
                        if NET1 is not None and NET2 is not None:
                            NET2_ratio_exp = (NET1 / NET2)**2
                            ratio_diff = abs(N0_ratio_meas - NET2_ratio_exp) / NET2_ratio_exp * 100
                            match = "YES" if ratio_diff < 10 else f"NO ({ratio_diff:.1f}%)"
                            
                            pair_str = f"{ch1_name} / {ch2_name}"
                            f.write(f"{pair_str:<50} {tel_yrs:<12.2f} {N0_ratio_meas:<18.4f} {NET2_ratio_exp:<18.4f} {match:<10}\n")
        
        # N_0 ratio test: same frequency, different tel_yrs
        f.write("\n\n" + "=" * 80 + "\n")
        f.write("N_0 RATIO TEST: Same frequency, different tel_yrs\n")
        f.write("Expected: N_0 ratio = (tel_yrs_2 / tel_yrs_1)\n")
        f.write("-" * 80 + "\n")
        f.write(f"{'Pair':<50} {'Freq (GHz)':<12} {'N0 ratio (meas)':<18} {'1/tel_yrs ratio (exp)':<22} {'Match?':<10}\n")
        f.write("-" * 112 + "\n")
        
        for freq, channels in sorted(freq_groups.items()):
            if len(channels) >= 2:
                channels_sorted = sorted(channels, key=lambda x: x[1])
                for i in range(len(channels_sorted)):
                    for j in range(i+1, len(channels_sorted)):
                        ch1_name, tel_yrs_1, N0_1 = channels_sorted[i]
                        ch2_name, tel_yrs_2, N0_2 = channels_sorted[j]
                        
                        N0_ratio_meas = N0_1 / N0_2
                        inv_telyrs_ratio_exp = tel_yrs_2 / tel_yrs_1
                        
                        ratio_diff = abs(N0_ratio_meas - inv_telyrs_ratio_exp) / inv_telyrs_ratio_exp * 100
                        match = "YES" if ratio_diff < 10 else f"NO ({ratio_diff:.1f}%)"
                        
                        pair_str = f"{ch1_name} / {ch2_name}"
                        f.write(f"{pair_str:<50} {freq:<12.0f} {N0_ratio_meas:<18.4f} {inv_telyrs_ratio_exp:<22.4f} {match:<10}\n")
        
        f.write("\n" + "=" * 80 + "\n")
    
    print(f"\nResults saved to: {results_file}")
    
    print(f"\n{'='*60}")
    print(f"Validation complete! Figures saved to: {validation_dir}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
