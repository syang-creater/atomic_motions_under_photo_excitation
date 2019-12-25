import sys;
from systemstructure import SystemStructure;
from displacement import Displacement;
import numpy as np;

class NormalCoordinate(Displacement):
    # /Users/shanyang/Desktop/SeSn-0.025-pnma/eigenvector.txt
    def read_eigenvector(self,atomic_mass,supercell_size,wave_vector,path):
        self._supercell_size = supercell_size;
        self._wave_vector = wave_vector;
        self._atomic_mass = dict(zip(atomic_mass,self._element_numbers))
        self._eigenvector_fract = np.loadtxt(path)
        self._eigenvector_direct_mass = np.empty((self._total_number,3),dtype=complex)
        count_start = 0;
        count_current = 0;
        for key in self._atomic_mass:
            value = self._atomic_mass[key]
            count_start = count_current;
            count_current = count_current+int(value);
            for atom in range(count_start, count_current):
                self._eigenvector_direct_mass[atom]=np.sqrt(self._total_number*np.prod(self._supercell_size))*float(key)*self._eigenvector_fract[atom]*\
                [np.linalg.norm(self._a_supercell_vector), np.linalg.norm(self._b_supercell_vector), np.linalg.norm(self._c_supercell_vector)]*\
                np.exp(2*np.pi*np.array([-1j])*np.dot(self._wave_vector,self._xdatcar_fract_unitcell[atom]))
        return self._eigenvector_direct_mass;

    def normalcoordinate(self):
        self._normal_coordinate_item = np.empty((self._count),dtype=complex)
        self._normal_coordinate = np.empty((self._steps),dtype=complex)
        for t in range(self._steps):
            for atom in range(self._total_number):
                self._normal_coordinate_item[t*self._total_number+atom] = np.dot(self._eigenvector_direct_mass[atom,:],self._displacement_direct_t[t*self._total_number+atom,:])\
                *np.exp(2*np.pi*np.array([-1j])*np.dot(self._wave_vector,self._xdatcar_fract_unitcell[t*self._total_number+atom]))
            self._normal_coordinate[t] = np.sum(self._normal_coordinate_item[t*self._total_number:(t+1)*self._total_number])/np.sqrt(np.prod(self._supercell_size)*self._total_number)
        return self._normal_coordinate

    def normalcoordinate_FFT(self,timestep_size,time_start, time_end):
        self._fs_nc = timestep_size;
        self._freq_nc = np.fft.fftfreq(time_end-time_start,self._fs_nc)
        self._normal_coordinate_fft  = np.abs(np.fft.fft(self._normal_coordinate[time_start:time_end]))
        return self._freq_nc, self._normal_coordinate_fft
