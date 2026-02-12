#!/usr/bin/env python
"""
Generate noise realizations for all SO channels specified in a YAML config file.
"""

import os
import yaml
import numpy as np
import healpy as hp
from tqdm import tqdm
from generate_noise_lat import SimonsObservatoryNoise

# Configuration
NSIMS = 10
NSPLITS = None
YAML_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/instr_params_lat_channels.yaml'
BASE_OUTPUT_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/output'

# Noise generation method: 'harmonic' or 'variance_map'
NOISE_METHOD = 'variance_map'
SURVEY = 'wide'  # Only used for variance_map method to determine which variance maps to use

# Directory containing variance maps (used if NOISE_METHOD='variance_map')
VARIANCE_MAP_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/variance_maps'
sohits_file = '/pscratch/sd/s/shamikg/so_lat_mapbased_noise/resources/so_mf_relhits_2048.fits'
sofoot_file = '/pscratch/sd/s/shamikg/so_lat_mapbased_noise/resources/so_common_lf_mf_uhf_bin_2048.fits'


def get_output_folder(yaml_path, base_output_dir):
    """Extract output folder name from last two underscore-separated words of yaml filename."""
    yaml_basename = os.path.basename(yaml_path)
    yaml_name = os.path.splitext(yaml_basename)[0]  # Remove .yaml extension
    parts = yaml_name.split('_')
    if len(parts) >= 2:
        folder_name = '_'.join(parts[-2:])
    else:
        folder_name = yaml_name
    return os.path.join(base_output_dir, folder_name)


def main():
    # Load channel names from YAML
    with open(YAML_FILE, 'r') as f:
        params = yaml.safe_load(f)
    
    channels = list(params.keys())
    print(f"Found {len(channels)} channels: {channels}")
    print(f"Noise generation method: {NOISE_METHOD}")
    if NOISE_METHOD == 'variance_map':
        print(f"Variance map directory: {VARIANCE_MAP_DIR}")

    # Create output folder
    output_folder = get_output_folder(YAML_FILE, BASE_OUTPUT_DIR)
    os.makedirs(output_folder, exist_ok=True)
    print(f"Output folder: {output_folder}")

    # Load hits and footprint maps for harmonic method
    so_hits = None
    so_foot = None
    if NOISE_METHOD == 'harmonic':
        print(f"Loading hits map: {sohits_file}")
        so_hits = hp.read_map(sohits_file)
        print(f"Loading footprint map: {sofoot_file}")
        so_foot = hp.read_map(sofoot_file)

    # Generate noise for each channel
    for ch_name in channels:
        print(f"\nProcessing channel: {ch_name}")
        
        # Initialize noise generator for this channel
        noise_gen = SimonsObservatoryNoise(
            ch_name, 
            params=YAML_FILE,
            noise_method=NOISE_METHOD,
            survey=SURVEY,
            so_hits=so_hits,
            so_foot=so_foot,
            variance_map_dir=VARIANCE_MAP_DIR
        )
        nside = noise_gen.nside

        for sim_idx in tqdm(range(NSIMS), desc=f"  {ch_name}", ncols=120):
            # Generate noise realization
            if NSPLITS is None:
                noise_map = noise_gen.get_noise()
                
                header = [
                    ('UNITS', 'uK_CMB', 'Map units'),
                    ('CHANNEL', noise_gen.channel, 'Channel name'),
                    ('FREQ', noise_gen.freq, 'Frequency in GHz'),
                    ('ELLKNEE_I', noise_gen.ell_knee_i, 'Intensity Knee multipole'),
                    ('ALPHA_I', noise_gen.alpha_i, 'Intensity Knee slope'),
                    ('ELLKNEE_P', noise_gen.ell_knee_p, 'Pol Knee multipole'),
                    ('ALPHA_P', noise_gen.alpha_p, 'Pol Knee slope'),
                    ('SIMIDX', sim_idx, 'Simulation index'),
                    ('NMETHOD', noise_gen.noise_method, 'Noise generation method'),
                ]
                
                # Add noise level if using harmonic method
                if noise_gen.noise_method == 'harmonic':
                    header.append(('NOISE', noise_gen.uKarcmin, 'noise level in uK-arcmin'))
                
                # Add telescope-years if using variance_map method
                if noise_gen.noise_method == 'variance_map' and noise_gen.tel_yrs is not None:
                    header.append(('TEL_YRS', noise_gen.tel_yrs, 'Telescope-years'))

                # Create output filename
                outfile = os.path.join(
                    output_folder,
                    f"sobs_lat-noise_{ch_name}_mc{sim_idx:03d}_nside{nside:04d}.fits"
                )

                # Save to FITS
                hp.write_map(outfile, noise_map, extra_header=header, overwrite=True, dtype=np.float32)
            else: 
                for split in range(NSPLITS):
                    noise_map = noise_gen.get_noise() * NSPLITS
                
                    header = [
                        ('UNITS', 'uK_CMB', 'Map units'),
                        ('CHANNEL', noise_gen.channel, 'Channel name'),
                        ('FREQ', noise_gen.freq, 'Frequency in GHz'),
                        ('ELLKNEE_I', noise_gen.ell_knee_i, 'Intensity Knee multipole'),
                        ('ALPHA_I', noise_gen.alpha_i, 'Intensity Knee slope'),
                        ('ELLKNEE_P', noise_gen.ell_knee_p, 'Pol Knee multipole'),
                        ('ALPHA_P', noise_gen.alpha_p, 'Pol Knee slope'),
                        ('SIMIDX', sim_idx, 'Simulation index'),
                        ('SPLIT', split, 'Split index'),
                        ('NMETHOD', noise_gen.noise_method, 'Noise generation method'),
                    ]
                    
                    # Add noise level if using harmonic method
                    if noise_gen.noise_method == 'harmonic':
                        header.append(('NOISE', noise_gen.uKarcmin, 'noise level in uK-arcmin'))
                    
                    # Add telescope-years if using variance_map method
                    if noise_gen.noise_method == 'variance_map' and noise_gen.tel_yrs is not None:
                        header.append(('TEL_YRS', noise_gen.tel_yrs, 'Telescope-years'))

                    # Create output filename
                    outfile = os.path.join(
                        output_folder,
                        f"sobs_lat-noise_{ch_name}_mc{sim_idx:03d}_split{split:02d}_nside{nside:04d}.fits"
                    )

                    # Save to FITS
                    hp.write_map(outfile, noise_map, extra_header=header, overwrite=True, dtype=np.float32)
        del noise_gen  # Free memory
        
    print(f"\nDone! All noise maps saved to {output_folder}")


if __name__ == '__main__':
    main()
