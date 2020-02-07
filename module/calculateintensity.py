import systemstructure;
import intensity;
import numpy as np;
import pdb

signal = intensity.Intensity();
signal.read_poscar('./test_example/POSCAR');
signal.read_xdatcar('/Users/shanyang/Desktop/SeSn-0.025-pnma/AIMD/old-XDATCAR/XDATCAR');
signal.get_info()
signal.warrap_error_fract()
signal.fract_supercell_to_unitcell()

wave_vector=np.array([3,0,1])
atomic_form_factor= np.array([34,50])
supercell_size = np.array([2,4,4])
timestep_size = 0.001 #(ps)
time_end = 12022;
time_start = 0;
signal.intensity(wave_vector,supercell_size,atomic_form_factor)
pdb.set_trace()
signal.plot_intensity(np.array([time_start,4253]))
signal.intensity_FFT(timestep_size,time_start,4253)
signal.plot_intensity_FFT(np.array([0,5]),np.array([0,10**5]))

signal.plot_intensity(np.array([4253,time_end]))
signal.intensity_FFT(timestep_size,4253,time_end)
signal.plot_intensity_FFT(np.array([0,5]),np.array([0,10**5]))

signal.intensity_auto_correlation(wave_vector,supercell_size,atomic_form_factor,time_start,4253)
signal.plot_intensity_auto_correlation(np.array([time_start,4253]))
signal.intensity_auto_correlation_FFT(timestep_size,time_start,4253)
signal.plot_intensity_auto_correlation_FFT(np.array([0,5]),np.array([0,10**5]))

signal.intensity_auto_correlation(wave_vector,supercell_size,atomic_form_factor,4253,time_end)
signal.plot_intensity_auto_correlation(np.array([4253,time_end]))
signal.intensity_auto_correlation_FFT(timestep_size,4253,time_end)
signal.plot_intensity_auto_correlation_FFT(np.array([0,5]),np.array([0,10**5]))
