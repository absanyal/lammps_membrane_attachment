# Polymer Simulation

This project simulates a polymer system using the LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) software.

## Files

- `polymer.py`: Generates the necessary input files for LAMMPS, including the data file that describes the polymer system and the input file that specifies the simulation parameters.
- `link_dist.py`: Analyzes the results of the simulation. It calculates the distance of the linkers from the membrane, detects when linkers hit the membrane, calculates the degree of attachment of the linkers to the membrane, and counts the number of linkers attached to the membrane. It can also plot these data and write them to text files.
- `run.sh`: Runs the simulation and the analysis.
- `submit_pc.sh`: Submits multiple simulation jobs.

## Usage

1. Run `polymer.py` to generate the LAMMPS input files.
```sh
python3 polymer.py
```
2. Run the LAMMPS simulation with the generated input file.
```sh
mpirun -n 8 lmp_mpi -i input.lammps > out.run
```

3. Run `link_dist.py` to analyze the results of the simulation.

## Requirements

- Python 3
- LAMMPS
- numpy
- matplotlib

## Details
- `polymer.py`: The `polymer.py` script sets up a simulation of a polymer system. It defines the parameters of the system, such as the number of atoms, bonds, and angles, the types of atoms and bonds, the dimensions of the simulation box, and the parameters for the LAMMPS simulation. The script also sets up the positions of the atoms in the polymer chain and the linkers, and the bonds and angles between them. The positions, bonds, and angles are written to the LAMMPS data file. The LAMMPS input file is also written by the script. It includes the commands to run the simulation, including the settings for the simulation (such as the timestep and the number of steps), the potentials for the interactions between the atoms, and the commands to output the results of the simulation.

#### Note
Please adjust the parameters in the polymer.py script to suit your specific needs before running the simulation.