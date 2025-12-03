import numpy as np
import healpy as hp
import yaml

ss = np.random.SeedSequence(6789012543567)

# Default SO parameters
DEFAULT_PARAMS = {
    'LF027': {'central_freq_GHz': 27., 'beam_fwhm_arcmin': 91., 'noise_uK_arcmin': 16., 'ell_knee': 30, 'alpha_knee': -2.4, 'nside': 512},
    'LF039': {'central_freq_GHz': 39., 'beam_fwhm_arcmin': 63., 'noise_uK_arcmin': 10., 'ell_knee': 30, 'alpha_knee': -2.4, 'nside': 512},
    'MF093': {'central_freq_GHz': 93., 'beam_fwhm_arcmin': 30., 'noise_uK_arcmin': 1.7, 'ell_knee': 50, 'alpha_knee': -2.5, 'nside': 512},
    'MF145': {'central_freq_GHz': 145., 'beam_fwhm_arcmin': 17., 'noise_uK_arcmin': 2.1, 'ell_knee': 50, 'alpha_knee': -3., 'nside': 512},
    'HF225': {'central_freq_GHz': 225., 'beam_fwhm_arcmin': 11., 'noise_uK_arcmin': 5.9, 'ell_knee': 70, 'alpha_knee': -3., 'nside': 512},
    'HF280': {'central_freq_GHz': 280., 'beam_fwhm_arcmin': 9., 'noise_uK_arcmin': 15., 'ell_knee': 100, 'alpha_knee': -3., 'nside': 512},
}

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

def _load_params(channel, params=None):
    """Load channel parameters from dict, yaml file path, or defaults."""
    if params is None:
        if channel not in DEFAULT_PARAMS:
            raise ValueError(f"Channel {channel} not found in default params. Must be one of {so_channels}")
        return DEFAULT_PARAMS[channel]
    
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
    def __init__(self, channel, params=None) -> None:
        self.channel = channel
        channel_params = _load_params(channel, params)
        
        self.uKarcmin = channel_params['noise_uK_arcmin']
        self.freq = channel_params['central_freq_GHz']
        self.ell_knee = channel_params['ell_knee']
        self.alpha = channel_params['alpha_knee']
        self.nside = channel_params['nside']
        self.npix = hp.nside2npix(self.nside)
        self.lmax = 3*self.nside - 1

        ALM = hp.Alm()
        self.almsz = ALM.getsize(self.lmax)
        self.zerom_idx = ALM.getidx(self.lmax, np.arange(self.lmax+1, dtype=int), m=0)

        ells = np.arange(self.lmax+1)
        atm_rescale = np.zeros(ells.shape)
        atm_rescale[2:] = 1. + (ells[2:] / self.ell_knee)**self.alpha
        
        self.Nl = uKarcmin2Nl(self.uKarcmin) * atm_rescale

        so_hits = hp.read_map(sohits_file)
        self.so_foot = hp.read_map(sofoot_file)

        self.hits_scaling = np.zeros(so_hits.shape)
        self.hits_scaling[so_hits>=1e-2] = 1. / np.sqrt(so_hits[so_hits>=1e-2])

        stream_idx = so_channels.index(channel)
        self.rng = np.random.default_rng(child_seeds[stream_idx])

    def get_noise(self, nside_out=None):
        hsqrt = np.sqrt(2.) / 2. 
        nlm_TEB = np.zeros((3, self.almsz), dtype=np.complex128)
        nlm_TEB[:, self.zerom_idx] = self.rng.normal(size=(3,len(self.zerom_idx)))
        nlm_TEB[:, self.zerom_idx[-1]:] = hsqrt * (self.rng.normal(size=(3,(self.almsz - len(self.zerom_idx)+1))) + 1j*self.rng.normal(size=(3,(self.almsz - len(self.zerom_idx)+1))))

        nlm_TEB[0] = hp.almxfl(nlm_TEB[0], np.sqrt(self.Nl)) / np.sqrt(2)
        nlm_TEB[1] = hp.almxfl(nlm_TEB[1], np.sqrt(self.Nl))
        nlm_TEB[2] = hp.almxfl(nlm_TEB[2], np.sqrt(self.Nl))

        noise_IQU = hp.alm2map(nlm_TEB, self.nside, pol=True) * self.so_foot * self.hits_scaling 

        return noise_IQU