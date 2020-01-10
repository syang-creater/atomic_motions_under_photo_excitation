import numpy as np;
from systemstructure import SystemStructure;

class Eigenvector(SystemStructure):
    def read_eigenvector(self,path):
        """
            the eigenvector belongs to the yaml file in phonopy, which contains the real part and imiginary part
        """
        self._read_eigenvector_unitcell=[];
        try:
            with open(path,'r') as eigenvector:
                while True:
                    line=eigenvector.readline();
                    if not line:
                        break;
                    else:
                        self._read_eigenvector_unitcell.append(line.strip().split());
        except FileNotFoundError as e:
            print('Eigenvector file does not exist: {}'.format(e))
        self._eigenvector_unitcell= np.asarray(self._read_eigenvector_unitcell);
        return;

    def generate_eigenvector_supercell(self,supercell_size):
        """
            transfer the unitcell eigenvector to supercell eigenvector
        """
        self._supercell_size = supercell_size;
        self._atoms_unitcell = self._eigenvector_unitcell.shape[0];
        self._eigenvector_supercell = [];
        for atom in range(self._atoms_unitcell):
            for x in range(self._supercell_size[0]):
                for y in range(self._supercell_size[1]):
                    for z in range(self._supercell_size[2]):
                        self._eigenvector_supercell.append(self._eigenvector_unitcell[atom]);
        self._eigenvector_supercell = np.asarray(self._eigenvector_supercell,dtype= float);
        return self._eigenvector_supercell;

    # def read_sposcar(self,supercell_path):
    #     self._sposcar = super().read_poscar(supercell_path)
    #     return self._sposcar;

    def eigenvector_supercell_phase(self,wave_vector):
        """
            the eigenvector in supercell does not include the relative phase of atoms in difference cell, which depends on the wave vector
            the wave vector is based on the unitcell
            self._coef:
            1. calculate the dot product to change the supercell fract to unitcell fract
            2. calculate the product of q*r
            3. calculate the phase exp(i*2*pi*q*r)
            comment: q is wave vector
        """
        # PATH:/Users/shanyang/Desktop/SeSn-0.025-pnma/infile.ssposcar
        self._wave_vector = wave_vector;
        self._coef=np.empty([self._eigenvector_supercell.shape[0],1],dtype=complex)
        self._eigenvector_supercell_phase=np.empty([self._eigenvector_supercell.shape[0],3],dtype=complex)
        for atom in range(self._total_number):
            self._coef[atom]= np.exp(2j*np.pi*np.dot(self._sposcar_fract[atom,:]*self._supercell_size,self._wave_vector))
            self._complex_eigen_x = np.complex(self._eigenvector_supercell[atom,0],self._eigenvector_supercell[atom,1]);
            self._complex_eigen_y = np.complex(self._eigenvector_supercell[atom,2],self._eigenvector_supercell[atom,3]);
            self._complex_eigen_z = np.complex(self._eigenvector_supercell[atom,4],self._eigenvector_supercell[atom,5]);
            self._eigenvector_supercell_phase[atom,:] = self._coef[atom]*np.array([self._complex_eigen_x,self._complex_eigen_y,self._complex_eigen_z],dtype=complex);
        #print(self._eigenvector_supercell_phase);
        return self._eigenvector_supercell_phase;

    def conj_eigenvector_supercell_phase_mass(self):
        """
            calculate all the parts conj(e)*np.sqrt(mass)
        """
        self._conj_eigenvector_supercell_phase_mass=np.empty([self._eigenvector_supercell.shape[0],3],dtype=complex)
        count_start = 0;
        count_current = 0;
        for key,value in zip(self._element_mass,self._element_numbers):
            count_start = count_current;
            count_current = count_current+int(value);
            for atom in range(count_start, count_current):
                self._conj_eigenvector_supercell_phase_mass[atom,:]= np.sqrt(key)*np.conj(self._eigenvector_supercell_phase[atom,:])
        return self._conj_eigenvector_supercell_phase_mass;
