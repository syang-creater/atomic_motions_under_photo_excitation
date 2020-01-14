import normalcoordinate;
import numpy as np;
import sys;
import matplotlib.pyplot as plt;
import pdb;

def main():
    """
        Calcuate normal coordinate for wave vector q and eigenvector at q for branch v
        the atomic order follow the TDEP for SPOSCAR and XDATCAR
    """
    nc=normalcoordinate.NormalCoordinate();
    nc.read_poscar('./test_example/POSCAR')
    nc.read_xdatcar('/Users/shanyang/Desktop/SeSn-0.025-pnma/AIMD/old-XDATCAR/XDATCAR')
    nc.get_info()
    nc.warrap_error_fract()

    nc.supercell_fract_to_direct()
    nc.fract_supercell_to_unitcell()
    nc.displacement_fract_at_t()
    nc.displacement_direct_at_t()

    supercell_size=np.array([2,4,4])
    atomic_mass= np.array([78.96,118.71])
    wave_vector=np.array([0,0,0])

    # load eigenvector
    # 0.3688768784THz TA at [0.25,0,0]
    nc.read_eigenvector('./test_example/eigenvector_TO_000')
    # load sposcar
    nc.read_sposcar('./test_example/infile.ssposcar');

    nc.generate_eigenvector_supercell(supercell_size)
    nc.eigenvector_supercell_phase(wave_vector)
    nc.atomic_mass(atomic_mass[0],atomic_mass[1])
    nc.conj_eigenvector_supercell_phase_mass()
    nc.normalcoordinate()
    nc.normalcoordinate_auto_correlate()
    time_start =0;
    time_end = 12022;

    #pdb.set_trace(); clear, step
    nc.plot_normalcoordinate(time_start, time_end)
    nc.plot_normalcoordinate_auto_correlate(time_start, time_end)

    timestep_size = 0.001
    #before photo_exciation
    nc.normalcoordinate_FFT(timestep_size, time_start, 4253)
    nc.normalcoordinate_auto_correlate_FFT(timestep_size, time_start, 4253)
    # xlim in THz
    nc.plot_normalcoordinate_FFT(np.array([0,5]));
    nc.plot_normalcoordinate_auto_correlate_FFT(np.array([0,5]));

    #after photo_excitation
    nc.normalcoordinate_FFT(timestep_size, 4253, time_end)
    nc.normalcoordinate_auto_correlate_FFT(timestep_size, 4253, time_end)

    # xlim in THz
    nc.plot_normalcoordinate_FFT(np.array([0,5]));
    nc.plot_normalcoordinate_auto_correlate_FFT(np.array([0,5]));


if __name__ == '__main__':
    main()
