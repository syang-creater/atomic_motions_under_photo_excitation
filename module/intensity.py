import sys;
from systemstructure import SystemStructure;
import matplotlib.pyplot as plt;
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

    def intensity_auto_correlation(self, wave_vector,supercell_size, atomic_form_factor, time_start, time_end):
        """
            t_0 is the time origin
            t_d is the time delay
            C_(t_d)=sum A(t_0)A(t_0+t_d) over from 0 to t_0 and take the average
        """
        self._intensity_auto_correlation = np.empty(time_end-time_start, dtype = complex)
        temp_intensity= self._structure_factor[time_start:time_end]
        for t in range(time_end-time_start):
            self._intensity_auto_correlation[t]= np.dot(temp_intensity,np.roll(temp_intensity,t))/(time_end-time_start);
        return self._intensity_auto_correlation.real

    def intensity_auto_correlation_FFT(self,timestep_size,time_start, time_end):
        """ calculate the FFT of intensity calculated from auto correlaion with timestep_size in unit of (ps), the unit of output frequency is (THz)"""
        self._fs = timestep_size;
        self._freq = np.fft.fftfreq(time_end-time_start,self._fs)
        self._intensity_auto_correlation_fft  = np.abs(np.fft.fft(self._intensity_auto_correlation.real))
        return self._freq, self._intensity_auto_correlation_fft

    def plot_intensity(self,xlim):
        plt.figure()
        plt.plot(self._intensity)
        plt.xlim(xlim[0],xlim[1])
        plt.xlabel('Time (fs)')
        plt.ylabel('Intensity')
        plt.show()
        return

    def plot_intensity_auto_correlation(self,xlim):
        plt.figure()
        plt.plot(np.arange(xlim[0],xlim[1]),self._intensity_auto_correlation.real)
        plt.xlim(xlim[0],xlim[1])
        plt.xlabel('Time (fs)')
        plt.ylabel('Intensity from auto correlation function')
        plt.show()
        return

    def plot_intensity_FFT(self,xlim,ylim):
        plt.figure()
        plt.plot(self._freq,self._fft)
        plt.xlim(xlim[0],xlim[1])
        plt.ylim(ylim[0],ylim[1])
        plt.xlabel('Frequency (THz)')
        plt.ylabel('Amplitude')
        plt.show()
        return

    def plot_intensity_auto_correlation_FFT(self,xlim,ylim):
        plt.figure()
        plt.plot(self._freq,self._intensity_auto_correlation_fft)
        plt.xlim(xlim[0],xlim[1])
        plt.ylim(ylim[0],ylim[1])
        plt.xlabel('Frequency (THz)')
        plt.ylabel('Amplitude')
        plt.show()
        return
