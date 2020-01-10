#!/Users/shanyang/Applications/anaconda3/bin/python
import systemstructure;
import displacement;
import intensity;
import normalcoordinate;
import numpy as np;

# all the input files follows the same atomic order and they are following the same atomic order as in TDEP

def main():
    """ import POSCAR SPOSCAR and XDARCAR """
    #AIMD=systemstructure.SystemStructure();
    #AIMD.read_poscar('/Users/shanyang/Desktop/SnS/harmonic-thermal-expansition/energy-cutoff/SnS-Pnma/POSCAR')
    #AIMD.read_xdatcar('/Users/shanyang/Desktop/SnS/harmonic-thermal-expansition/energy-cutoff/SnS-Pnma/MD-cutoff/MD-1/XDATCAR')
    """ output system informations """
    #AIMD.get_info();
    #AIMD.atomic_mass(10,15);

    #AIMD.warrap_error_fract()
    #AIMD.supercell_fract_to_direct()

    #AIMD.fract_supercell_to_unitcell()
    #AIMD.warrap_error_fract_unitcell()

    """ Calculate atomic displacement """
    #disp=displacement.Displacement();
    #disp.read_poscar('/Users/shanyang/Desktop/SnS/harmonic-thermal-expansition/energy-cutoff/SnS-Pnma/POSCAR')
    #disp.read_xdatcar('/Users/shanyang/Desktop/SnS/harmonic-thermal-expansition/energy-cutoff/SnS-Pnma/MD-cutoff/MD-1/XDATCAR')
    #disp.get_info()
    #disp.warrap_error_fract()

    #disp.supercell_fract_to_direct()
    #disp.displacement_fract_at_t()
    #disp.displacement_direct_at_t()

    #disp.fract_supercell_to_unitcell()
    #disp.warrap_error_fract_unitcell()
    #disp.warrap_error_fract_unitcell()

    """ Calculate Bragg intensity """
    #i = intensity.Intensity();
    #i.read_poscar('/Users/shanyang/Desktop/SnS/harmonic-thermal-expansition/energy-cutoff/SnS-Pnma/POSCAR');
    #i.read_xdatcar('/Users/shanyang/Desktop/SnS/harmonic-thermal-expansition/energy-cutoff/SnS-Pnma/MD-cutoff/MD-1/XDATCAR');
    #i.get_info()
    #i.warrap_error_fract()
    #i.fract_supercell_to_unitcell()

    #wave_vector=np.array([2,4,3])
    #atomic_form_factor= np.array([30,50])
    #supercell_size = np.array([2,4,4])
    #timestep_size = 0.001 #(ps)
    #i.intensity(wave_vector,supercell_size,atomic_form_factor)
    #print(i.intensity_FFT(timestep_size,0,255))

    """ Calcuate normal coordinate """
    nc=normalcoordinate.NormalCoordinate();
    nc.read_poscar('/Users/shanyang/Desktop/SnS/harmonic-thermal-expansition/energy-cutoff/SnS-Pnma/POSCAR')
    nc.read_xdatcar('/Users/shanyang/Desktop/SnS/harmonic-thermal-expansition/energy-cutoff/SnS-Pnma/MD-cutoff/MD-1/XDATCAR')
    nc.get_info()
    nc.warrap_error_fract()

    nc.supercell_fract_to_direct()
    nc.fract_supercell_to_unitcell()
    # in order to calculate direct coordinate, need to calculate fract coordinate first
    nc.displacement_fract_at_t()
    nc.displacement_direct_at_t()

    # load eigenvector
    supercell_size=np.array([2,4,4])
    atomic_mass= np.array([45,80])
    wave_vector=np.array([0,0,0]) # the wave vector that you used to generate eigenvector
    nc.read_eigenvector('/Users/shanyang/Desktop/photo_excitation_module/module/eigenvector_example')
    nc.read_sposcar('/Users/shanyang/Desktop/SeSn-0.025-pnma/infile.ssposcar');
    nc.generate_eigenvector_supercell(supercell_size)
    nc.eigenvector_supercell_phase(wave_vector)
    nc.atomic_mass(atomic_mass[0],atomic_mass[1])
    nc.conj_eigenvector_supercell_phase_mass()
    nc.normalcoordinate()

    timestep_size = 0.001 #(ps)
    nc.normalcoordinate_FFT(timestep_size, 0, 255)


if __name__ == '__main__':
    main()
