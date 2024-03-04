import numpy as np
import angle
from angle import filter_angle

############ ENTER INPUTS HERE #######################

# various toggles
auto_generate_seed = 1
dump_minimization = 1

# Name of polymer data file
data_fname_str = 'single.data'

# Name of LAMMPS input file
input_fname_str = 'single_input.lammps'

# Name of information file
# info_fname_str = 'data/info.txt'

# Name of linker position file
single_pos_fname_str = 'data/single_pos.txt'


# ---Types---
atom_types = 1

# ---LAMMPS INPUT FILE PARAMETERS START---
# These directly get written into the LAMMPS file

units = 'lj'
dimension = 3
boundary = 'p p p'

atom_style = 'molecular'

global_cutoff = 0.0

timestep = 0.00001

run_steps = 5000000

measure_distance_every = 100

# Brownian parameters
brn_T = 310
brn_gamma = 1.0
brn_seed = 490563

if (auto_generate_seed == 1):
    brn_seed = np.random.randint(1, 100000000)

print("Seed used =", brn_seed)

brownian = [brn_T, brn_seed, brn_gamma]

# ---LAMMPS INPUT FILE PARAMETERS END---

n_atoms = 1

# ---Box dimensions---
xlo, xhi = 0.0, 1000
ylo, yhi = 0.0, 1000
zlo, zhi = 0.0, 1000

# ---Setup mass---
mass = [
    [1, 1.0, "monomer_chain1"],
]

################## END INPUTS #####################

# ---Setup positions---
positions = []

chain = 1
normalatom = 1
for i in range(n_atoms):
    thisatom = normalatom

    # For chain parallel to surface
    px = (xhi - xlo)/2.0
    py = (yhi - ylo)/2.0
    pz = (zhi - zlo)/2.0

    positions.append([chain, thisatom, px, py, pz])


############## DATA FILES ################

# ---Write data file for atoms---
with open(data_fname_str, 'w') as data_f:

    # ---Header---
    data_f.write("Two chains and floating linkers\n\n")

    # Numbers
    data_f.write('{} atoms\n'.format(n_atoms))

    data_f.write('\n')

    # Types

    data_f.write('{} atom types\n'.format(atom_types))

    data_f.write('\n')

    # Box size
    data_f.write('{} {} xlo xhi\n'.format(xlo, xhi))
    data_f.write('{} {} ylo yhi\n'.format(ylo, yhi))
    data_f.write('{} {} zlo zhi\n'.format(zlo, zhi))

    data_f.write('\n')

    # Masses
    data_f.write('Masses \n\n')

    for i in range(atom_types):
        data_f.write(
            '{} {} # {}\n'.format(mass[i][0], mass[i][1], mass[i][2]))

    data_f.write('\n')

    # ---Atoms---
    data_f.write('Atoms\n\n')

    for i, pos in enumerate(positions):
        data_f.write('{} {} {} {} {} {}\n'.format(i+1, *pos))

    data_f.write('\n')

# ---Write LAMMPS input file---
with open(input_fname_str, 'w') as input_f:

    input_f.write('units {}\n'.format(units))
    input_f.write('dimension {}\n'.format(dimension))
    input_f.write('boundary {}\n\n'.format(boundary))

    input_f.write('atom_style {}\n\n'.format(atom_style))

    input_f.write('read_data {}\n\n'.format(data_fname_str))

    input_f.write('timestep {0:.10f}\n\n'.format(timestep))
    
    input_f.write('pair_style zero {} nocoeff\n'.format(global_cutoff))
    input_f.write('pair_coeff * *\n\n')

    input_f.write('variable x equal x[1]\n')
    input_f.write('variable y equal y[1]\n')
    input_f.write('variable z equal z[1]\n\n')

    input_f.write('variable tsteps equal time\n\n')

    input_f.write(
        'fix link_pos all print 100 "${tsteps} ${x} ${y} ${z}"')
    input_f.write(' file {} screen no\n\n'.format(single_pos_fname_str))

    input_f.write('dump mydump all atom 1000 dump.single.lammpstrj\n\n')

    input_f.write(
        'fix 1 all brownian {0} {1} gamma_t {2}\n\n'.format(*brownian))

    input_f.write('thermo_style custom step time temp etotal\n')
    input_f.write('thermo 100000\n\n')

    input_f.write('run {}\n'.format(run_steps))
