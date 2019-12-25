import sys;
from systemstructure import SystemStructure;
import numpy as np;

class Intensity(SystemStructure):
    def intensity(self, wave_vector,supercell_size, atomic_form_factor):
        """ calculate bragg intensity for every step for wave_vector based on unitcell """
        self._element_form_factor = dict(zip(atomic_form_factor,self._element_numbers))
        self._structure_factor = np.empty(self._steps,dtype=complex);
        self._intensity = np.empty(self._steps, dtype = float);
        temp_matrix = np.empty(self._count,dtype=complex);
        for t in range(self._steps):
            count_start = 0;
            count_current = 0;
            for key in self._element_form_factor:
                value = self._element_form_factor[key]
                count_start = count_current;
                count_current = count_current+int(value);
                for atom in range(count_start, count_current):
                    temp_matrix[t*self._total_number+atom] = float(key)*np.exp(2*np.pi*np.array([-1j])*np.dot(self._xdatcar_fract_unitcell[t*self._total_number+atom,:],wave_vector))
            self._structure_factor[t] = np.sum(temp_matrix[t*self._total_number:(t+1)*self._total_number], axis=0)/np.prod(supercell_size)
            self._intensity[t] = np.absolute(self._structure_factor[t])**2
        return self._intensity

    def intensity_FFT(self,timestep_size,time_start, time_end):
        """ calculate the FFT of intensity with timestep_size in unit of (ps), the unit of output frequency is (THz)"""
        self._fs = timestep_size;
        self._freq = np.fft.fftfreq(time_end-time_start,self._fs)
        self._fft  = np.abs(np.fft.fft(self._intensity[time_start:time_end]))
        return self._freq, self._fft
