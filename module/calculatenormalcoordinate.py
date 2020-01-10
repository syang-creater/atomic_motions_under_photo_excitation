import normalcoordinate;
import numpy as np;
import sys;
import matplotlib.pyplot as plt;

def main():
    """
        Calcuate normal coordinate for wave vector q and eigenvector at q for branch v
        the atomic order follow the TDEP for SPOSCAR and XDATCAR
    """
    nc=normalcoordinate.NormalCoordinate();
    nc.read_poscar('./test_example/POSCAR')
    nc.read_xdatcar('./test_example/XDATCAR-315')
    nc.get_info()
    nc.warrap_error_fract()

    nc.supercell_fract_to_direct()
    nc.fract_supercell_to_unitcell()
    nc.displacement_fract_at_t()
    nc.displacement_direct_at_t()

    supercell_size=np.array([2,4,4])
    atomic_mass= np.array([78.96,118.71])
    wave_vector=np.array([0.25,0,0])

    # load eigenvector
    # 0.3688768784THz TA at [0.25,0,0]
    nc.read_eigenvector('./test_example/eigenvector_TA_02500')
    # load sposcar
    nc.read_sposcar('./test_example/infile.ssposcar');

    nc.generate_eigenvector_supercell(supercell_size)
    nc.eigenvector_supercell_phase(wave_vector)
    nc.atomic_mass(atomic_mass[0],atomic_mass[1])
    nc.conj_eigenvector_supercell_phase_mass()
    nc.normalcoordinate()
    nc.normalcoordinate_square()
    time_start =0;
    time_end = 8500;
    nc.plot_normalcoordinate(time_start, time_end)
    nc.plot_normalcoordinate_square(time_start, time_end)

    timestep_size = 0.001
    nc.normalcoordinate_FFT(timestep_size, time_start, time_end-6500)
    nc.normalcoordinate_square_FFT(timestep_size, time_start, time_end-6500)
    # xlim in THz
    nc.plot_normalcoordinate_FFT(np.array([0,5]));
    nc.plot_normalcoordinate_square_FFT(np.array([0,5]));

    #timestep_size = 0.001
    nc.normalcoordinate_FFT(timestep_size, 2000, time_end)
    nc.normalcoordinate_square_FFT(timestep_size, 2000, time_end)
    # xlim in THz
    nc.plot_normalcoordinate_FFT(np.array([0,5]));
    nc.plot_normalcoordinate_square_FFT(np.array([0,5]));


if __name__ == '__main__':
    main()
