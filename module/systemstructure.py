import sys;
import numpy as np;

class SystemStructure:
    def read_poscar(self,path):
        """ read the POSCAR (VASP) file, which you used to generate supercell """
        try:
            with open(path,'r') as poscar:
                self._system_unitcell=poscar.readline()
                self._scale_unitcell=float(poscar.readline().rstrip('\n').rstrip());
                self._a_unitcell_vector=np.array([float(i)*self._scale_unitcell for i in poscar.readline().rstrip('\n').split()])
                self._b_unitcell_vector=np.array([float(i)*self._scale_unitcell for i in poscar.readline().rstrip('\n').split()])
                self._c_unitcell_vector=np.array([float(i)*self._scale_unitcell for i in poscar.readline().rstrip('\n').split()])
                self._latticevector_matrix_unitcell=np.round(np.stack((self._a_unitcell_vector,self._b_unitcell_vector,self._c_unitcell_vector)),6)
                poscar.readline()
                poscar.readline()
                self._poscar=[]
                while True:
                    line=poscar.readline().rstrip('\n').split();
                    if not line:
                        break
                    if (self._isfloat(*[items for items in line])):
                        self._poscar.append(line)
                self._poscar_fract = np.asarray(self._poscar, dtype= float)
        except FileNotFoundError as e:
            print('POSCAR file does not exist: {}'.format(e))
        return;

    def read_sposcar(self,path):
        """ read the SPOSCAR (VASP) suprcell file """
        try:
            with open(path,'r') as sposcar:
                sposcar.readline()
                sposcar.readline()
                sposcar.readline()
                sposcar.readline()
                sposcar.readline()
                self._element_names = [name for name in sposcar.readline().rstrip('\n').split()]
                self._element_numbers = np.array([int(number) for number in sposcar.readline().rstrip('\n').split()])
                self._total_number = np.sum(self._element_numbers)
                self._sposcar=[]
                while True:
                    line=sposcar.readline().rstrip('\n').split();
                    if not line:
                        break
                    if (self._isfloat(*[items for items in line])):
                        self._sposcar.append(line)
                self._sposcar_fract = np.asarray(self._sposcar, dtype= float)
        except FileNotFoundError as e:
            print('SPOSCAR file does not exist: {}'.format(e))
        except IOE as e:
            print('IO error'.format(e))
            raise e
        return self._sposcar_fract;

    def read_xdatcar(self,path):
        """ read the XDARCAR file from AIMD (VASP) """
        try:
            with open(path,'r') as xdatcar:
                self._system=xdatcar.readline()
                self._scale_supercell=float(xdatcar.readline().rstrip('\n').rstrip());
                self._a_supercell_vector=np.array([float(i)*self._scale_supercell for i in xdatcar.readline().rstrip('\n').split()])
                self._b_supercell_vector=np.array([float(i)*self._scale_supercell for i in xdatcar.readline().rstrip('\n').split()])
                self._c_supercell_vector=np.array([float(i)*self._scale_supercell for i in xdatcar.readline().rstrip('\n').split()])
                self._latticevector_matrix_supercell=np.round(np.stack((self._a_supercell_vector,self._b_supercell_vector,self._c_supercell_vector)),6)
                self._element_names = [name for name in xdatcar.readline().rstrip('\n').split()]
                self._element_numbers = np.array([int(number) for number in xdatcar.readline().rstrip('\n').split()])
                self._total_number = np.sum(self._element_numbers)
                self._xdatcar=[]
                self._count = 0
                while True:
                    line=xdatcar.readline().rstrip('\n').split();
                    if not line:
                        break
                    if (self._isfloat(*[items for items in line])):
                        self._xdatcar.append(line)
                        self._count +=1
                self._xdatcar_fract = np.asarray(self._xdatcar,dtype = float)
                self._steps = int(self._count/self._total_number)
        except FileNotFoundError as e:
            print('XDARCAR file does not exist:{}'.format(e))
            raise e
        except IOE as e:
            print('IO error'.format(e))
            raise e
        return self._xdatcar_fract;

    def _isfloat(self,*value):
        for it in value:
            try:
                float(it)
            except ValueError:
                return False
        return True;

    def get_info(self):
        """ Print out system informations from POSCAR and XDATCAR file."""
        for self._name, self._number in zip(self._element_names,self._element_numbers):
            print('Element name: %s Element number: %s' %(self._name,self._number))
        print('Total number: %s' %(self._total_number))
        print('Time steps: %s' %(self._steps))
        print('Lattice parameter for unitcell:')
        print(self._latticevector_matrix_unitcell)
        print('Lattice parameter for supercell:')
        print(self._latticevector_matrix_supercell)
        return

    def atomic_mass(self,*arg):
        """take atomic mass in element order from input unit (a.u) """
        self._element_mass=[]
        for mass in arg:
            self._element_mass.append(float(mass));
        self._element_mass = np.asarray(self._element_mass, dtype = float)

    def warrap_error_fract(self):
        """ clean warrap error when atoms are jumping from one side to the other """
        self._equil_fract = self._xdatcar_fract[0:self._total_number,:]
        for t in range(int(self._steps)):
            for atom in range(self._total_number):
                for xyz in range(3):
                    if (self._xdatcar_fract[t*self._total_number+atom,xyz]-self._equil_fract[atom,xyz]) > 0.5:
                        self._xdatcar_fract[t*self._total_number+atom,xyz]=self._xdatcar_fract[t*self._total_number+atom,xyz]-1
                    elif (self._xdatcar_fract[t*self._total_number+atom,xyz]-self._equil_fract[atom,xyz]) < (-0.5):
                        self._xdatcar_fract[t*self._total_number+atom,xyz]=self._xdatcar_fract[t*self._total_number+atom,xyz]+1
        return self._xdatcar_fract

    def supercell_fract_to_direct(self):
        """ charge from fract coordinate to direct coordinate unit(A)"""
        self._xdatcar_direct_supercell = np.zeros((self._count,3))
        abc_supercell=[np.linalg.norm(self._a_supercell_vector), np.linalg.norm(self._b_supercell_vector), np.linalg.norm(self._c_supercell_vector)]
        for posi in range(self._count):
            self._xdatcar_direct_supercell[posi,:] = self._xdatcar_fract[posi,:]*abc_supercell
        return self._xdatcar_direct_supercell

    def fract_supercell_to_unitcell(self):
        """ change fraction coordinate from supercell to unitcell """
        self._xdatcar_fract_unitcell = np.zeros((self._count,3))
        abc_supercell=[np.linalg.norm(self._a_supercell_vector), np.linalg.norm(self._b_supercell_vector), np.linalg.norm(self._c_supercell_vector)]
        abc_unitcell=[np.linalg.norm(self._a_unitcell_vector), np.linalg.norm(self._b_unitcell_vector), np.linalg.norm(self._c_unitcell_vector)]
        for posi in range(self._count):
            self._xdatcar_fract_unitcell[posi,:] = self._xdatcar_fract[posi,:]*abc_supercell/abc_unitcell%np.array([1,1,1])
        return self._xdatcar_fract_unitcell

    def warrap_error_fract_unitcell(self):
        """ clean warrap error when atoms are jumping from one unitcell to the other unitcell """
        self._equil_fract_unitcell = self._xdatcar_fract_unitcell[0:self._total_number,:]
        for t in range(int(self._steps)):
            for atom in range(self._total_number):
                for xyz in range(3):
                    if (self._xdatcar_fract_unitcell[t*self._total_number+atom,xyz]-self._equil_fract_unitcell[atom,xyz]) >= 0.5:
                        self._xdatcar_fract_unitcell[t*self._total_number+atom,xyz]=self._xdatcar_fract_unitcell[t*self._total_number+atom,xyz]-1
                    elif (self._xdatcar_fract_unitcell[t*self._total_number+atom,xyz]-self._equil_fract_unitcell[atom,xyz]) <= (-0.5):
                        self._xdatcar_fract_unitcell[t*self._total_number+atom,xyz]=self._xdatcar_fract_unitcell[t*self._total_number+atom,xyz]+1
        return self._xdatcar_fract_unitcell
