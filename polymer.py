import numpy as np

############ ENTER INPUTS HERE #######################

auto_generate_seed = 1

# Name of polymer data file
data_fname_str = 'polymer.data'

# Name of LAMMPS input file
input_fname_str = 'input.lammps'

# Name of information file
info_fname_str = 'info.txt'

# ---Types---
atom_types = 2
bond_types = 1
angle_types = 1

# ---LAMMPS INPUT FILE PARAMETERS START---
# These directly get written into the LAMMPS file

units = 'lj'
dimension = 3
boundary = 'p p p'

atom_style = 'molecular'

global_cutoff = 10.0

timestep = 0.00001

thermalize_steps = 10000
run_steps = 3000000

measure_distance_every = 10000

# Brownian parameters
brn_T = 310
brn_gamma = 1.0
brn_seed = 490563

if (auto_generate_seed == 1):
    brn_seed = np.random.randint(1, 100000000)

print("Seed used =", brn_seed)

brownian = [brn_T, brn_seed, brn_gamma]

# Potentials in format epsilon, sigma, Rc

wall_chain = [5.0, 5.0, 5.0]
wall_linker = [1500.0, 2.5, 10.0]

# Format: Bond_style name, bond type number, k, r0
bonds_styles = [
    ['harmonic', 1, 1500.0, 2.5]
]

# Format: angle_style name, angle type number, k, theta0
angle_styles = [
    ['harmonic', 1, 15000.0, 177.55]
]

# ---LAMMPS INPUT FILE PARAMETERS END---


# ---Filament Parameters

# Distance between two atoms in the filament
bondlength = 2.5

# Chain info (only count polymer chain)
n_chains = 1
chain_offset = 10
distance_from_axis = 500

# Per chain numbers
n_atoms = 20
n_bonds = n_atoms - 1
n_angles = n_bonds - 1

n_cross_bonds = 0

# ---Linker numbers---
n_skip_mem_linkers = 4  # Gap between two successive linkers
n_linkers_cross = 0

# ---Box dimensions---
xlo, xhi = 0.0, 1000
ylo, yhi = 0.0, 350
zlo, zhi = 0.0, 350

# ---Setup mass---
mass = [
    [1, 1.0, "monomer_chain1"],
    [2, 1.1, "linker_membrane"]
]

################## END INPUTS #####################

# ---Setup linker numbers---
n_linkers_membrane = 0
for i in range(n_atoms):
    if (i % n_skip_mem_linkers == 0):
        n_linkers_membrane += 1

n_bonds += n_linkers_membrane

# ---Setup positions---
positions = []

chain = 1
normalatom = 1
for i in range(n_atoms):
    thisatom = normalatom

    px = (xhi - xlo)/2.0
    # px = 100
    py = (yhi - ylo)/2.0 + distance_from_axis
    pz = (i * bondlength) - (zhi - zlo)/2
    positions.append([chain, thisatom, px, py, pz])

# # Linkers
thisatom = 2
chain = 1
print("membrane linkers =", n_linkers_membrane)
for i in range(n_linkers_membrane):
    j = i * n_skip_mem_linkers
    # print(i, j)
    px = positions[j][2]
    py = positions[j][3] + 3.0
    pz = positions[j][4]
    positions.append([chain, thisatom, px, py, pz])


# ---Setup bonds----
bonds = []

# linear bonds in chain1, bond type = 1
bond_type = 1
for i in range(n_atoms - 1):
    b_start = i+1
    b_stop = b_start + 1
    bond = [bond_type, b_start, b_stop]
    bonds.append(bond)

for i in range(n_linkers_membrane):
    b_start = i + 1 + n_atoms
    b_stop = i * n_skip_mem_linkers + 1
    bond = [bond_type, b_start, b_stop]
    bonds.append(bond)

# ---Setup angles---
angles = []

# 180 degree angle between two successive bonds, chain 1, angle type=1
angle_type = 1
for i in range(n_atoms - 2):
    a_1 = i+1
    a_2 = a_1 + 1
    a_3 = a_2 + 1
    angle = [angle_type, a_1, a_2, a_3]
    angles.append(angle)
    

############## DATA FILES ################

# ---Write data file for information---
with open(info_fname_str, 'w') as info_f:
    info_f.write('{}\n'.format(n_atoms))
    info_f.write('{}\n'.format(n_linkers_membrane))

# ---Write data file for atoms---
with open(data_fname_str, 'w') as data_f:

    # ---Header---
    data_f.write("Two chains and floating linkers\n\n")

    # Numbers
    data_f.write('{} atoms\n'.format(n_atoms * n_chains +
                                     n_linkers_membrane + n_linkers_cross))
    data_f.write('{} bonds\n'.format(n_bonds * n_chains + n_cross_bonds))
    data_f.write('{} angles\n'.format(n_angles * n_chains))

    data_f.write('\n')

    # Types

    data_f.write('{} atom types\n'.format(atom_types))
    data_f.write('{} bond types\n'.format(bond_types))
    data_f.write('{} angle types\n'.format(angle_types))

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

    # ---Bonds---
    data_f.write('Bonds\n\n')

    for i, bond in enumerate(bonds):
        data_f.write('{} {} {} {}\n'.format(i+1, *bond))

    data_f.write('\n')

    # ---Angles---
    data_f.write('Angles\n\n')

    for i, angle in enumerate(angles):
        data_f.write('{} {} {} {} {}\n'.format(i+1, *angle))

# ---Write LAMMPS input file---
with open(input_fname_str, 'w') as input_f:

    input_f.write('units {}\n'.format(units))
    input_f.write('dimension {}\n'.format(dimension))
    input_f.write('boundary {}\n\n'.format(boundary))

    input_f.write('atom_style {}\n\n'.format(atom_style))

    input_f.write('read_data {}\n\n'.format(data_fname_str))

    input_f.write('region membrane cylinder x {} {} {} 0 {}\n\n'.format(
        yhi/2, zhi/2, yhi/2, xhi))

    input_f.write('group chain1 type 1\n'.format())
    input_f.write('group link_mem type 2\n\n'.format())

    input_f.write('pair_style zero {} nocoeff\n'.format(global_cutoff))
    input_f.write('pair_coeff * *\n\n')

    for bst in bonds_styles:
        input_f.write('bond_style {}\n'.format(bst[0]))
        input_f.write('bond_coeff {} {} {}\n\n'.format(*bst[1:]))

    for ast in angle_styles:
        input_f.write('angle_style {}\n'.format(ast[0]))
        input_f.write('angle_coeff {} {} {}\n\n'.format(*ast[1:]))

    input_f.write('timestep {0:.10f}\n\n'.format(timestep))

    input_f.write('dump minimization all atom 1 dump.min.lammpstrj\n')
    input_f.write('minimize 0.0 1.0e-8 10000 10000\n\n')

    input_f.write('fix 1 all nve/limit 0.01\n\n')

    input_f.write(
        'fix wallchain1 chain1 wall/region membrane lj93 {} {} {}\n'.format(*wall_chain))
    input_f.write(
        'fix walllinkermem link_mem wall/region membrane lj93 {} {} {}\n\n'.format(*wall_linker))

    for i in range(n_linkers_membrane):
        j = i + n_atoms + 1

        input_f.write('variable x{} equal x[{}]\n'.format(i+1, j))
        input_f.write('variable y{} equal y[{}]\n'.format(i+1, j))
        input_f.write('variable z{} equal z[{}]\n\n'.format(i+1, j))

    input_f.write('thermo_style custom step time temp etotal\n')
    input_f.write('thermo 10000\n\n')

    input_f.write('run {}\n\n'.format(10000))

    input_f.write('unfix 1\n')
    input_f.write('undump minimization\n\n')
    
    input_f.write('reset_timestep 0\n\n')
    
    input_f.write('variable tsteps equal time\n\n')

    input_f.write('fix link_pos all print {} "${{tsteps}} '.format(measure_distance_every))

    for i in range(n_linkers_membrane):
        input_f.write('${{x{0}}} ${{y{0}}} ${{z{0}}} '.format(i+1))

    input_f.write('" file link_pos.txt screen no\n\n')

    input_f.write('dump mydump all atom 1000 dump.lammpstrj\n\n')

    input_f.write(
        'fix 2 all brownian {0} {1} gamma_t {2}\n\n'.format(*brownian))

    input_f.write('thermo 100000\n\n')

    input_f.write('run {}\n'.format(run_steps))