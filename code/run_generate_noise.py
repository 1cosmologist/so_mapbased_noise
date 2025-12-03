#!/usr/bin/env python
"""
Generate noise realizations for all SO channels specified in a YAML config file.
"""

import os
import yaml
import numpy as np
import healpy as hp
from tqdm import tqdm
from generate_noise import SimonsObservatoryNoise

# Configuration
NSIMS = 5
NSPLITS = None
YAML_FILE = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/instr_params_baseline_pessimistic.yaml'
BASE_OUTPUT_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/output'


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

    # Create output folder
    output_folder = get_output_folder(YAML_FILE, BASE_OUTPUT_DIR)
    os.makedirs(output_folder, exist_ok=True)
    print(f"Output folder: {output_folder}")

    # Generate noise for each channel
    for ch_name in channels:
        print(f"\nProcessing channel: {ch_name}")
        
        # Initialize noise generator for this channel
        noise_gen = SimonsObservatoryNoise(ch_name, params=YAML_FILE)
        nside = noise_gen.nside

        for sim_idx in tqdm(range(NSIMS), desc=f"  {ch_name}", ncols=120):
            # Generate noise realization
            if NSPLITS is None:
                noise_map = noise_gen.get_noise()
                
                header = [
                    ('UNITS', 'uK_CMB', 'Map units'),
                    ('CHANNEL', noise_gen.channel, 'Channel name'),
                    ('FREQ', noise_gen.freq, 'Frequency in GHz'),
                    ('NOISE', noise_gen.uKarcmin, 'noise level in uK-arcmin'),
                    ('ELLKNEE', noise_gen.ell_knee, 'Knee multipole'),
                    ('ALPHA', noise_gen.alpha, 'Knee slope'),
                    ('SIMIDX', sim_idx, 'Simulation index'),
                ]

                # Create output filename
                outfile = os.path.join(
                    output_folder,
                    f"sobs_noise_{ch_name}_mc{sim_idx:03d}_nside{nside:04d}.fits"
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
                        ('NOISE', noise_gen.uKarcmin, 'noise level in uK-arcmin'),
                        ('ELLKNEE', noise_gen.ell_knee, 'Knee multipole'),
                        ('ALPHA', noise_gen.alpha, 'Knee slope'),
                        ('SIMIDX', sim_idx, 'Simulation index'),
                        ('SPLIT', split, 'Split index'),
                    ]

                    # Create output filename
                    outfile = os.path.join(
                        output_folder,
                        f"sobs_noise_{ch_name}_mc{sim_idx:03d}_split{split:02d}_nside{nside:04d}.fits"
                    )

                    # Save to FITS
                    hp.write_map(outfile, noise_map, extra_header=header, overwrite=True, dtype=np.float32)
        del noise_gen  # Free memory
        
    print(f"\nDone! All noise maps saved to {output_folder}")


if __name__ == '__main__':
    main()
