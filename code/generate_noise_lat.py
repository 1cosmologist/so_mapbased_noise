import numpy as np
import healpy as hp
import yaml
import os
import glob

ss = np.random.SeedSequence(1234567890123)

# Default SO parameters
DEFAULT_PARAMS = {
    'LF027': {'central_freq_GHz':  27., 'beam_fwhm_arcmin': 7.4, 'noise_uK_arcmin': 44., 'ell_knee_i': 500,  'alpha_knee_i': -3.5, 'ell_knee_p': 700, 'alpha_knee_p': -1.4, 'nside': 2048},
    'LF039': {'central_freq_GHz':  39., 'beam_fwhm_arcmin': 5.1, 'noise_uK_arcmin': 23., 'ell_knee_i': 500,  'alpha_knee_i': -3.5, 'ell_knee_p': 700, 'alpha_knee_p': -1.4, 'nside': 2048},
    'MF093': {'central_freq_GHz':  93., 'beam_fwhm_arcmin': 2.2, 'noise_uK_arcmin': 3.8, 'ell_knee_i': 2100, 'alpha_knee_i': -3.5, 'ell_knee_p': 700, 'alpha_knee_p': -1.4, 'nside': 2048},
    'MF145': {'central_freq_GHz': 145., 'beam_fwhm_arcmin': 1.4, 'noise_uK_arcmin': 4.1, 'ell_knee_i': 3000, 'alpha_knee_i': -3.5, 'ell_knee_p': 700, 'alpha_knee_p': -1.4, 'nside': 2048},
    'HF225': {'central_freq_GHz': 225., 'beam_fwhm_arcmin': 1.0, 'noise_uK_arcmin': 10., 'ell_knee_i': 3800, 'alpha_knee_i': -3.5, 'ell_knee_p': 700, 'alpha_knee_p': -1.4, 'nside': 2048},
    'HF280': {'central_freq_GHz': 280., 'beam_fwhm_arcmin': 0.9, 'noise_uK_arcmin': 25., 'ell_knee_i': 3800, 'alpha_knee_i': -3.5, 'ell_knee_p': 700, 'alpha_knee_p': -1.4, 'nside': 2048},
}
# Channel name to frequency mapping for variance map file lookup
CHANNEL_FREQ_MAP = {
    'LF027': 'f027', 'LF039': 'f039',
    'MF093': 'f093', 'MF145': 'f145',
    'HF225': 'f225', 'HF280': 'f280', 
}

# Default variance maps directory
DEFAULT_VARIANCE_MAP_DIR = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/variance_maps'

so_channels = list(DEFAULT_PARAMS.keys())

nfreqs_tot = 25 #len(so_channels)  # Budget for 25 frequency channels
# Generate deterministic integer seeds from the SeedSequence for each channel
# This ensures reproducibility and avoids SeedSequence state consumption issues
_child_ss = ss.spawn(nfreqs_tot)
child_seeds = [cs.generate_state(1)[0] for cs in _child_ss]

sohits_file = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/so_sat_relhits_C_nside512.fits'
sofoot_file = '/pscratch/sd/s/shamikg/so_mapbased_noise/resources/so_sat_full-binary_C_nside512.fits'

def uKarcmin2Nl(uKarcmin):
    return (uKarcmin * np.deg2rad(1./60.))**2

def find_variance_map(channel, tel_yrs=None, variance_map_dir=None, survey='wide'):
    """
    Find variance map file for a given channel.
    
    Parses files named like: so_f027_noise_var_2.50tel-yrs.fits
    
    Parameters
    ----------
    channel : str
        Channel name (e.g., 'LF027', 'MF093')
    tel_yrs : float, optional
        Number of telescope-years for this channel. If provided, looks for exact match.
        If None, returns the first available variance map for the channel.
    variance_map_dir : str, optional
        Directory to search for variance maps
        
    Returns
    -------
    str or None
        Path to the variance map if found, None otherwise.
    """
    if variance_map_dir is None:
        variance_map_dir = DEFAULT_VARIANCE_MAP_DIR
    
    freq_code = CHANNEL_FREQ_MAP.get(channel)
    if freq_code is None:
        raise ValueError(f"Unknown channel {channel}. Must be one of {list(CHANNEL_FREQ_MAP.keys())}")
    
    # Look for exact match with specified tel-yrs
    if survey == 'wide': filename = f"so_lat_wide_1_el_{freq_code}_noise_var_{tel_yrs:.2f}tel-yrs.fits"
    if survey == 'delens-wide': filename = f"so_lat_delens-wide_1_el_{freq_code}_noise_var_{tel_yrs:.2f}tel-yrs.fits"
    if survey == 'delens-bk': filename = f"so_lat_delens-bk_1_el_{freq_code}_noise_var_{tel_yrs:.2f}tel-yrs.fits"
    
    filepath = os.path.join(variance_map_dir, filename)
    if os.path.exists(filepath):
        return filepath
    else:
        raise FileNotFoundError(f"No variance map found for channel {channel} with tel_yrs={tel_yrs:.2f} at {filepath}")

def _load_params(channel, params=None):
    """Load channel parameters from dict, yaml file path, or defaults."""
    if params is None:
        if channel not in DEFAULT_PARAMS:
            raise ValueError(f"Channel {channel[:5]} not found in default params. Must be one of {so_channels}")
        return DEFAULT_PARAMS[channel[:5]]
    
    if isinstance(params, str):
        # Assume it's a yaml file path
        with open(params, 'r') as f:
            params = yaml.safe_load(f)
    
    if isinstance(params, dict):
        if channel in params:
            return params[channel]
        else:
            raise ValueError(f"Channel {channel} not found in provided params")
    
    raise TypeError("params must be a dict or a path to a yaml file")
    
class SimonsObservatoryNoise:
    """
    Generate noise realizations for Simons Observatory channels.
    
    Parameters
    ----------
    channel : str
        Channel name (e.g., 'LF027', 'MF093', etc.)
    params : dict or str, optional
        Channel parameters as dict or path to YAML file. For variance_map method,
        the config should include 'tel_yrs' for each channel.
    noise_method : str, optional
        Method for noise generation:
        - 'harmonic': Generate noise in harmonic space with hits-based inhomogeneity (default)
        - 'variance_map': Generate noise using per-pixel variance maps
    so_hits : ndarray, optional
        Relative hits map. Required for harmonic method.
    so_foot : ndarray, optional
        Binary footprint map. Required for harmonic method.
    variance_map_path : str, optional
        Path to variance map file. If None and noise_method='variance_map',
        will auto-detect from DEFAULT_VARIANCE_MAP_DIR using tel_yrs from config
    variance_map_dir : str, optional
        Directory to search for variance maps (used if variance_map_path is None)
    """
    def __init__(self, channel, params=None, noise_method='harmonic', survey='wide',
                 so_hits=None, so_foot=None,
                 variance_map_path=None, variance_map_dir=None) -> None:
        self.channel = channel
        self.noise_method = noise_method
        self.survey = survey
        channel_params = _load_params(channel, params)
        
        # Get tel_yrs from config if present (required for variance_map method)
        self.tel_yrs = channel_params.get('tel_yrs', None)
        
        # noise_uK_arcmin is only required for harmonic method
        self.uKarcmin = channel_params.get('noise_uK_arcmin', None)
        self.freq = channel_params['central_freq_GHz']
        
        # Parse params for I and P
        self.ell_knee_i = channel_params.get('ell_knee_i', 3000)
        self.alpha_i = channel_params.get('alpha_knee_i', -3.5)
        self.ell_knee_p = channel_params.get('ell_knee_p', 700)
        self.alpha_p = channel_params.get('alpha_knee_p', -1.4)
        
        self.ell_knee = self.ell_knee_p  # Legacy support
        self.alpha = self.alpha_p # Legacy support
        
        self.nside = channel_params['nside']
        self.npix = hp.nside2npix(self.nside)
        self.lmax = 3*self.nside - 1

        ALM = hp.Alm()
        self.almsz = ALM.getsize(self.lmax)
        self.zerom_idx = ALM.getidx(self.lmax, np.arange(self.lmax+1, dtype=int), m=0)

        ells = np.arange(self.lmax+1)
        
        # Precompute scalings for I and P
        self.ell_scaling_i = np.ones(ells.shape)
        self.ell_scaling_i[2:] = np.sqrt(1. + (ells[2:] / self.ell_knee_i)**self.alpha_i)
        
        self.ell_scaling_p = np.ones(ells.shape)
        self.ell_scaling_p[2:] = np.sqrt(1. + (ells[2:] / self.ell_knee_p)**self.alpha_p)
        
        # Legacy atmospheric rescaling (using P params)
        atm_rescale = np.zeros(ells.shape)
        atm_rescale[2:] = 1. + (ells[2:] / self.ell_knee)**self.alpha
        
        # Nl is only needed for harmonic method
        if noise_method == 'harmonic':
            if self.uKarcmin is None:
                raise ValueError(f"noise_uK_arcmin is required in config for channel {channel} when using harmonic method")
            self.Nl = uKarcmin2Nl(self.uKarcmin) * atm_rescale
        else:
            self.Nl = None
        
        # Ell scaling for 1/f noise (used in variance_map method)
        self.ell_scaling = np.zeros(ells.shape)
        self.ell_scaling[2:] = np.sqrt(atm_rescale[2:])

        # hits_scaling and so_foot only needed for harmonic method
        self.so_foot = None
        self.hits_scaling = None
        if noise_method == 'harmonic':
            if so_hits is None or so_foot is None:
                raise ValueError("so_hits and so_foot must be provided for harmonic noise method")
            self.so_foot = so_foot
            self.hits_scaling = np.zeros(so_hits.shape)
            self.hits_scaling[so_hits >= 1e-2] = 1. / np.sqrt(so_hits[so_hits >= 1e-2])
        
        # Load variance map if using variance_map method
        self.noise_variance_map = None
        self.variance_map_path = None
        if noise_method == 'variance_map':
            if variance_map_path is not None:
                self.variance_map_path = variance_map_path
            else:
                self.variance_map_path = find_variance_map(channel[:5], tel_yrs=self.tel_yrs, variance_map_dir=variance_map_dir, survey=self.survey)
            
            if self.variance_map_path is None:
                tel_yrs_msg = f" with tel_yrs={self.tel_yrs:.2f}" if self.tel_yrs is not None else ""
                raise FileNotFoundError(
                    f"No variance map found for channel {channel[:5]}{tel_yrs_msg}. "
                    f"Please provide variance_map_path or ensure files exist in {variance_map_dir or DEFAULT_VARIANCE_MAP_DIR}"
                )
            
            self.noise_variance_map = hp.read_map(self.variance_map_path)
            print(f"  Loaded variance map: {self.variance_map_path}")

        if self.survey == 'wide':
            stream_idx = so_channels.index(channel[:5])
        elif self.survey == 'delens-wide':
            stream_idx = so_channels.index(channel[:5]) + len(so_channels)
        elif self.survey == 'delens-bk':
            stream_idx = so_channels.index(channel[:5]) + 2*len(so_channels)
            
        self.rng = np.random.default_rng(child_seeds[stream_idx])

    def get_noise(self, nside_out=None):
        """
        Generate a noise realization.
        
        Uses the method specified during initialization:
        - 'harmonic': Generate alms with proper Nl, then apply hits-based inhomogeneity
        - 'variance_map': Generate white noise in pixel space, apply 1/f in harmonic space,
                          then scale by sqrt(variance_map)
        
        Parameters
        ----------
        nside_out : int, optional
            Output nside (not currently used, reserved for future)
            
        Returns
        -------
        noise_IQU : ndarray
            Shape (3, npix) noise map in I, Q, U
        """
        if self.noise_method == 'harmonic':
            return self._get_noise_harmonic()
        elif self.noise_method == 'variance_map':
            return self._get_noise_variance_map()
        else:
            raise ValueError(f"Unknown noise_method: {self.noise_method}. Must be 'harmonic' or 'variance_map'")
    
    def _get_noise_harmonic(self):
        """
        Generate noise using harmonic space method with hits-based inhomogeneity.
        
        This is the original method that generates alms with the proper noise power spectrum
        and applies spatial inhomogeneity via relative hits map.
        """
        hsqrt = np.sqrt(2.) / 2. 
        nlm_TEB = np.zeros((3, self.almsz), dtype=np.complex128)
        nlm_TEB[:, self.zerom_idx] = self.rng.normal(size=(3,len(self.zerom_idx)))
        nlm_TEB[:, self.zerom_idx[-1]:] = hsqrt * (self.rng.normal(size=(3,(self.almsz - len(self.zerom_idx)+1))) + 1j*self.rng.normal(size=(3,(self.almsz - len(self.zerom_idx)+1))))

        nlm_TEB[0] = hp.almxfl(nlm_TEB[0], np.sqrt(self.Nl)) / np.sqrt(2)
        nlm_TEB[1] = hp.almxfl(nlm_TEB[1], np.sqrt(self.Nl))
        nlm_TEB[2] = hp.almxfl(nlm_TEB[2], np.sqrt(self.Nl))

        noise_IQU = hp.alm2map(nlm_TEB, self.nside, pol=True) * self.so_foot * self.hits_scaling 

        return noise_IQU
    
    def _get_noise_variance_map(self):
        """
        Generate noise using per-pixel variance map method.
        
        Similar to the approach in so_350_depthmap.ipynb:
        1. Generate unit Gaussian noise in pixel space
        2. Apply 1/f ell scaling in harmonic space
        3. Multiply by sqrt(variance_map) for spatial inhomogeneity
        
        The variance map contains the noise variance per pixel in uK^2-pixel units,
        computed from NET and observation time.
        """
        # Generate unit Gaussian noise in pixel space (I, Q, U)
        noise_pix = self.rng.normal(loc=0.0, scale=1.0, size=(3, self.npix))
        
        # Transform to harmonic space and apply 1/f ell scaling
        nlm = hp.map2alm(noise_pix, lmax=self.lmax, pol=True)
        
        # Apply intensity scaling to I component
        hp.almxfl(nlm[0], self.ell_scaling_i, inplace=True)
        
        # Apply polarization scaling to Q/U components
        hp.almxfl(nlm[1], self.ell_scaling_p, inplace=True)
        hp.almxfl(nlm[2], self.ell_scaling_p, inplace=True)
        
        # Transform back to pixel space and apply variance map scaling
        noise_pix = hp.alm2map(nlm, nside=self.nside, pol=True) * np.sqrt(self.noise_variance_map)
        
        # Rescale polarization noise level to temperature 
        noise_pix[0] /= np.sqrt(2.)
        
        return noise_pix