python3 single_particle.py
time mpirun -n 1 lmp_mpi -i single_input.lammps > out.run
python3 D_single.py