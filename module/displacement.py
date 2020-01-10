import sys;
import numpy as np;
from systemstructure import SystemStructure;
import matplotlib.pyplot as plt

class Displacement(SystemStructure):
    """ calculate atomic displacement in unit (A)"""

    def displacement_fract_at_t(self):
        self._displacement_fract_t=np.zeros((self._count,3))
        self._xdatcar_fract_supercell_equal = self._xdatcar_fract[0:self._total_number,:]
        for t in range(self._steps):
            self._displacement_fract_t[t*self._total_number:(t+1)*self._total_number,:]=self._xdatcar_fract[t*self._total_number:(t+1)*self._total_number,:]-self._xdatcar_fract_supercell_equal
        return self._displacement_fract_t

    def displacement_direct_at_t(self):
        self._displacement_direct_t=np.zeros((self._count,3))
        abc_supercell=[np.linalg.norm(self._a_supercell_vector), np.linalg.norm(self._b_supercell_vector), np.linalg.norm(self._c_supercell_vector)]
        for t in range(self._count):
            self._displacement_direct_t[t*self._total_number:(t+1)*self._total_number,:]=self._displacement_fract_t[t*self._total_number:(t+1)*self._total_number,:]*abc_supercell
        return self._displacement_direct_t
