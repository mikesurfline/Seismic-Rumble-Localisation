import numpy as np
from skimage.feature import structure_tensor

def enh(spectrogram, sigma=(4,250)):
    A1,A2,A3 = structure_tensor(spectrogram,sigma=sigma,mode='nearest')
    A3 = (A3 - np.min(A3)) / (np.max(A3) - np.min(A3))
    return spectrogram*A3

spect_seis = {}
spect_seis['fs'] = 200
spect_seis['NFFT'] = 400
spect_seis['noverlap'] = 198
spect_seis['nperseg'] = 199
spect_seis['frequency_resolution'] = 2
spect_seis['time_resolution'] = 0.005
spect_seis['fft_resolution'] = spect_seis['NFFT']
spect_seis['fft_stride'] = spect_seis['NFFT'] - spect_seis['noverlap']

spect_shum = {}
spect_shum['fs'] = 8000
spect_shum['NFFT'] = 2000
spect_shum['noverlap'] = 1970 #1975 for good resolution
spect_shum['nperseg'] = 1999 #1999 for good resolution
spect_shum['fft_resolution'] = spect_shum['NFFT']
spect_shum['fft_stride'] = spect_shum['NFFT'] - spect_shum['noverlap']
