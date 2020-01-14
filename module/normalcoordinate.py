import sys;
from systemstructure import SystemStructure;
from displacement import Displacement;
from eigenvector import Eigenvector;
import matplotlib.pyplot as plt;

import numpy as np;

class NormalCoordinate(Displacement,Eigenvector):
    # /Users/shanyang/Desktop/SeSn-0.025-pnma/eigenvector.txt
    """
        N in the normal coordinate is the number of unitcells in supercell
    """
    def normalcoordinate(self):
        self._normal_coordinate_item = np.empty((self._count),dtype=complex)
        self._normal_coordinate = np.empty((self._steps),dtype=complex)
        for t in range(self._steps):
            for atom in range(self._total_number):
                self._normal_coordinate_item[t*self._total_number+atom] = np.dot(self._conj_eigenvector_supercell_phase_mass[atom,:],self._displacement_direct_t[t*self._total_number+atom,:])
            self._normal_coordinate[t] = np.sum(self._normal_coordinate_item[t*self._total_number:(t+1)*self._total_number])/np.prod(self._supercell_size)
        return self._normal_coordinate.real

    def normalcoordinate_auto_correlate(self):
        self._normal_coordinate_auto_correlate=np.correlate(self._normal_coordinate,self._normal_coordinate,'same')
        return self._normal_coordinate_auto_correlate

    def normalcoordinate_FFT(self,timestep_size,time_start, time_end):
        self._fs_nc = timestep_size;
        self._freq_nc = np.fft.fftfreq(time_end-time_start,self._fs_nc)
        self._normal_coordinate_fft  = np.abs(np.fft.fft(self._normal_coordinate.real[time_start:time_end]))
        return self._freq_nc, self._normal_coordinate_fft

    def normalcoordinate_auto_correlate_FFT(self,timestep_size,time_start, time_end):
        self._normal_coordinate_auto_correlate_fft  = np.abs(np.fft.fft(self._normal_coordinate_auto_correlate[time_start:time_end]))
        return self._freq_nc, self._normal_coordinate_auto_correlate_fft

    def plot_normalcoordinate(self,time_start,time_end,step=1):
        plt.figure()
        plt.plot(np.arange(time_start,time_end),self._normal_coordinate.real[time_start:time_end])
        plt.xlim([time_start, time_end])
        plt.xlabel('Time (fs)')
        plt.ylabel(r'Normal Coordinate Q ($amu^{0.5}\AA$)')
        plt.show()
        return

    def plot_normalcoordinate_auto_correlate(self,time_start,time_end,step=1):
        plt.figure()
        plt.plot(np.arange(time_start,time_end),self._normal_coordinate_auto_correlate[time_start:time_end])
        plt.xlim([time_start, time_end])
        plt.xlabel('Time (fs)')
        plt.ylabel(r'Normal Coordinate $Q^{2}$ ($amu\AA^{2}$)')
        plt.show()
        return

    def plot_normalcoordinate_FFT(self,xlim):
        plt.figure()
        plt.plot(self._freq_nc,self._normal_coordinate_fft)
        plt.xlim(xlim[0],xlim[1])
        plt.xlabel('Frequency (THz)')
        plt.ylabel('Amplitude')
        plt.xlim()
        plt.show()
        return

    def plot_normalcoordinate_auto_correlate_FFT(self,xlim):
        plt.figure()
        plt.plot(self._freq_nc,self._normal_coordinate_auto_correlate_fft)
        plt.xlim(xlim[0],xlim[1])
        plt.xlabel('Frequency (THz)')
        plt.ylabel('Amplitude')
        plt.xlim()
        plt.show()
        return
