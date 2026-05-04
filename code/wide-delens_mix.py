#!/usr/bin/env python
"""
Mix wide and delensing LAT noise maps using per-pixel inverse variance weighting
for a given telescope-year split fraction.

For a fraction f in [0, 1]:
  - Wide survey receives    (1-f) of the total telescope-years per channel
  - Delensing survey receives  f  of the total telescope-years per channel
  - f = 0  →  output is purely wide-survey noise
  - f = 1  →  output is purely delensing-survey noise

The combination is a per-pixel ILC (inverse variance weighting).  Allocating
only a fraction of the budget to a survey inflates its pixel variance relative
to the stored variance map by 1/fraction, so the stored noise realization is
rescaled by 1/sqrt(fraction) before combining.

Derivation (all operations per pixel):
  var_w_new  = var_w  / (1-f)        ← fewer tel-yrs → more noise
  var_d_new  = var_d  /  f
  n_w_new    = n_w    / sqrt(1-f)    ← rescale stored realisation
  n_d_new    = n_d    / sqrt(f)

  ILC combination (minimises per-pixel variance):
      n_comb = (sqrt(1-f)*var_d*n_w + sqrt(f)*var_w*n_d) / ((1-f)*var_d + f*var_w)

Pixels observed by only one survey use that survey's rescaled noise alone.
"""

import os
import yaml
import numpy as np
import healpy as hp
from tqdm import tqdm

# ── Configuration ────────────────────────────────────────────────────────────

NSIMS = 10

# Fraction of total telescope-years assigned to the delensing survey.
DELENS_FRACTION = 0.9

WIDE_YAML   = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/instr_params_lat-wide_channels.yaml'
DELENS_YAML = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/instr_params_lat-delens_channels.yaml'

WIDE_NOISE_DIR   = '/pscratch/sd/s/shamikg/so_mapbased_noise/output/lat-wide_channels'
DELENS_NOISE_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/output/lat-delens_channels'
VARIANCE_MAP_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/variance_maps'
BASE_OUTPUT_DIR  = '/pscratch/sd/s/shamikg/so_mapbased_noise/output'

# First-5-char channel name → frequency code used in variance-map filenames
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
    return hp.read_map(fpath)


def noise_filepath(noise_dir, ch_name, sim_idx, nside):
    """Return path to a pre-saved noise-map FITS file."""
    return os.path.join(
        noise_dir,
        f"sobs_lat-noise_{ch_name}_mc{sim_idx:03d}_nside{nside:04d}.fits",
    )


# ── Core mixing routine ───────────────────────────────────────────────────────

def combine_noise_maps(n_wide, n_delens, var_wide, var_delens, f):
    """
    Combine wide and delensing noise maps for a telescope-year split fraction f.

    Parameters
    ----------
    n_wide, n_delens : ndarray, shape (3, npix)
        Noise realisations (I, Q, U) generated with the full (yaml) tel-yrs budget.
    var_wide, var_delens : ndarray, shape (npix,)
        Per-pixel noise variance maps (uK^2-pixel) matching the yaml tel-yrs.
    f : float
        Fraction of total tel-yrs assigned to delensing (0 ≤ f ≤ 1).

    Returns
    -------
    n_comb : ndarray, shape (3, npix)
    """
    a = 1.0 - f   # wide fraction
    b = f          # delensing fraction

    w_in = var_wide   > 0   # pixels inside the wide footprint
    d_in = var_delens > 0   # pixels inside the delensing footprint

    n_comb = np.zeros_like(n_wide)

    # ── Pixels only observed by the wide survey ──────────────────────────────
    # Noise amplitude scales as 1/sqrt(tel-yrs), so rescale by 1/sqrt(a).
    w_only = w_in & ~d_in
    if np.any(w_only):
        if a == 0:
            pass  # f=1: wide gets 0 tel-yrs → no contribution
        else:
            n_comb[:, w_only] = n_wide[:, w_only] / np.sqrt(a)

    # ── Pixels only observed by the delensing survey ─────────────────────────
    d_only = ~w_in & d_in
    if np.any(d_only):
        if b == 0:
            pass  # f=0: delensing gets 0 tel-yrs → no contribution
        else:
            n_comb[:, d_only] = n_delens[:, d_only] / np.sqrt(b)

    # ── Pixels observed by both surveys → ILC ────────────────────────────────
    both = w_in & d_in
    if np.any(both):
        if a == 0:   # f=1 → pure delensing
            n_comb[:, both] = n_delens[:, both]
        elif b == 0: # f=0 → pure wide
            n_comb[:, both] = n_wide[:, both]
        else:
            vw = var_wide[both]    # (M,)
            vd = var_delens[both]  # (M,)
            nw = n_wide[:, both]   # (3, M)
            nd = n_delens[:, both] # (3, M)
            # ILC formula derived in module docstring; (M,) arrays broadcast over (3, M)
            denom = a * vd + b * vw                                    # (M,)
            n_comb[:, both] = (np.sqrt(a) * vd * nw
                                + np.sqrt(b) * vw * nd) / denom        # (3, M)

    return n_comb


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    f = float(DELENS_FRACTION)
    assert 0.0 <= f <= 1.0, "DELENS_FRACTION must be in [0, 1]"

    print(f"Delensing fraction  f = {f:.4f}  "
          f"(wide: {1-f:.4f}, delensing: {f:.4f})")

    # Load YAML configs
    with open(WIDE_YAML)   as fh: wide_params   = yaml.safe_load(fh)
    with open(DELENS_YAML) as fh: delens_params  = yaml.safe_load(fh)

    channels = list(wide_params.keys())
    if set(channels) != set(delens_params.keys()):
        print("Warning: channel lists differ — using intersection.")
        channels = sorted(set(channels) & set(delens_params.keys()))
    print(f"Channels ({len(channels)}): {channels}")

    output_dir = os.path.join(BASE_OUTPUT_DIR, f'lat-mix_f{f:.2f}')
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output directory: {output_dir}\n")

    for ch_name in channels:
        print(f"Channel: {ch_name}")
        wp = wide_params[ch_name]
        dp = delens_params[ch_name]
        nside      = wp['nside']
        w_tel_yrs  = wp['tel_yrs']
        d_tel_yrs  = dp['tel_yrs']
        print(f"  Wide yaml tel-yrs: {w_tel_yrs}  |  Delens yaml tel-yrs: {d_tel_yrs}")

        # Pre-load variance maps (same for all sims of this channel)
        var_wide   = load_variance_map(ch_name, w_tel_yrs, 'wide',       VARIANCE_MAP_DIR)
        var_delens = load_variance_map(ch_name, d_tel_yrs, 'delens_wide', VARIANCE_MAP_DIR)

        for sim_idx in tqdm(range(NSIMS), desc=f"  {ch_name}", ncols=120):
            nw = hp.read_map(
                noise_filepath(WIDE_NOISE_DIR,   ch_name, sim_idx, nside),
                field=[0, 1, 2], 
            )
            nd = hp.read_map(
                noise_filepath(DELENS_NOISE_DIR, ch_name, sim_idx, nside),
                field=[0, 1, 2], 
            )

            n_comb = combine_noise_maps(nw, nd, var_wide, var_delens, f)

            header = [
                ('UNITS',   'uK_CMB',              'Map units'),
                ('CHANNEL', ch_name,               'Channel name'),
                ('FREQ',    wp['central_freq_GHz'],'Frequency in GHz'),
                ('SIMIDX',  sim_idx,               'Simulation index'),
                ('MIXFRAC', f,                     'Delensing tel-yrs fraction (0=wide, 1=delens)'),
                ('W_TYRS',  w_tel_yrs,             'Wide survey yaml tel-yrs (full budget)'),
                ('D_TYRS',  d_tel_yrs,             'Delens survey yaml tel-yrs (full budget)'),
            ]
            outfile = os.path.join(
                output_dir,
                f"sobs_lat-noise-mix_{ch_name}_f{f:.2f}_mc{sim_idx:03d}_nside{nside:04d}.fits",
            )
            hp.write_map(outfile, n_comb, extra_header=header,
                         overwrite=True, dtype=np.float32)

    print(f"\nDone. Mixed noise maps saved to {output_dir}")


if __name__ == '__main__':
    main()
