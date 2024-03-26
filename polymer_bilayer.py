import numpy as np
import angle
from angle import filter_angle

############ ENTER INPUTS HERE #######################

# various toggles
auto_generate_seed = 1
dump_minimization = 1
make_walls_big = 0
compute_CoM_MSD = 0
compute_CoM_pos = 1
linkers_on_filament = 0

# Name of polymer data file
data_fname_str = 'polymer_bilayer.data'

# Name of LAMMPS input file
input_fname_str = 'input_bilayer.lammps'

# Name of information file
info_fname_str = 'data/info.txt'

# Name of linker position file
link_pos_fname_str = 'data/link_pos.txt'

# Name of end positions file
end_pos_fname_str = 'data/e2e_pos.txt'

# Name of CoM MSD file
com_msd_fname_str = 'data/com_msd.txt'

# Name of CoM position file
com_pos_fname_str = 'data/com_pos.txt'


# ---Types---
atom_types = 3

bond_types = 1
angle_types = 1

# ---LAMMPS INPUT FILE PARAMETERS START---
# These directly get written into the LAMMPS file

units = 'lj'
dimension = 3
boundary = 'p p p'

atom_style = 'molecular'

global_cutoff = 8.0

timestep = 0.00001

thermalize_steps = 1000
run_steps = 3000000

measure_distance_every = 10000

# Brownian parameters
brn_T = 310
brn_gamma = 1.0
brn_seed = 490563

brn_gamma_embedded = 100.0
brn_seed_embedded = 490563

if (auto_generate_seed == 1):
    brn_seed = np.random.randint(1, 100000000)
    brn_seed_embedded = np.random.randint(1, 100000000)

print("Seed used (filament) =", brn_seed)
print("Seed used (embedded) =", brn_seed_embedded)

brownian = [brn_T, brn_seed, brn_gamma]
brownian_embedded = [brn_T, brn_seed_embedded, brn_gamma_embedded]

# Potentials in format epsilon, sigma, Rc

wall_chain = [5.0, 1.0, 1.0]
wall_linker = [1500.0, 2.5, 10.0]
wall_memlinker = [5.0, 1.0, 1.0]
memlinker_chain = [1500.0, 2.5, 8.0]

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

# Angle between the chain and the membrane (in degrees)
theta = 45
theta = filter_angle(theta)
print("theta in degrees =", theta)

theta = np.radians(theta)
print("theta in radians =", theta)

buffer = 5

# Chain info (only count polymer chain)
n_chains = 1
chain_offset = 10
distance_from_axis = 0
distance_from_axis = 170

distance_between_membranes = 10

# distance_from_axis -= (buffer + distance_between_membranes)

# Per chain numbers
n_atoms = 20
n_bonds = n_atoms - 1
n_angles = n_bonds - 1

n_cross_bonds = 0

# ---Linker numbers---
n_skip_fil_linkers = 4  # Gap between two successive linkers
n_linkers_cross = 0
n_linkers_embedded = 1000

# ---Box dimensions---
xlo, xhi = 0.0, 1000
ylo, yhi = 0.0, 380
zlo, zhi = 0.0, 380

if (make_walls_big == 1):
    xlo, xhi = 0.0, 1000
    ylo, yhi = 0.0, 500
    zlo, zhi = 0.0, 500

# ---Setup mass---
mass = [
    [1, 1.0, "monomer_chain1"],
    [2, 1.1, "linker_filament"],
    [3, 1.0, "linker_embedded"]
]

################## END INPUTS #####################

# ---Setup linker numbers---
n_linkers_filament = 0
if (linkers_on_filament == 1):
    for i in range(n_atoms):
        if (i % n_skip_fil_linkers == 0):
            n_linkers_filament += 1

n_bonds += n_linkers_filament

# ---Setup positions---
positions = []

chain = 1
normalatom = 1
for i in range(n_atoms):
    thisatom = normalatom

    # # For chain parallel to surface
    # px = (xhi - xlo)/2.0
    # py = (yhi - ylo)/2.0 + distance_from_axis
    # pz = -(i * bondlength) + (zhi - zlo)/2
    
    # For chain perpendicular to surface/at an angle
    px = (xhi - xlo)/2.0
    py = (yhi - ylo)/2.0 + distance_from_axis - i * bondlength * np.cos(theta)
    pz = (zhi - zlo)/2 - i * bondlength * np.sin(theta)
    
    
    positions.append([chain, thisatom, px, py, pz])

# Filament Linkers
thisatom = 2
chain = 1
print("linkers on chain =", n_linkers_filament)
for i in range(n_linkers_filament):
    j = i * n_skip_fil_linkers
    j = n_atoms - j - 1
    print("L {} placed with m {}".format(i+1+n_atoms, j+1))
    px = positions[j][2] - bondlength
    py = positions[j][3]
    pz = positions[j][4]
    positions.append([chain, thisatom, px, py, pz])

# Embedded linkers
print("linkers embedded =", n_linkers_embedded)


# ---Setup bonds----
bonds = []

# linear bonds in chain1, bond type = 1
bond_type = 1
for i in range(n_atoms - 1):
    b_start = i+1
    b_stop = b_start + 1
    bond = [bond_type, b_start, b_stop]
    bonds.append(bond)

if (linkers_on_filament == 1):
    for i in range(n_linkers_filament):
        b_start = i + 1 + n_atoms
        b_stop = i * n_skip_fil_linkers + 1
        b_stop = n_atoms - b_stop + 1
        bond = [bond_type, b_start, b_stop]
        print("L {} attached to m {}".format(b_start, b_stop))
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
    info_f.write('{}\n'.format(n_linkers_filament))
    info_f.write('{}\n'.format(n_skip_fil_linkers))

# ---Write data file for atoms---
with open(data_fname_str, 'w') as data_f:

    # ---Header---
    data_f.write("Two chains and floating linkers\n\n")

    # Numbers
    data_f.write('{} atoms\n'.format(n_atoms * n_chains +
                                     n_linkers_filament + n_linkers_cross))
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
    
    outer_radius = (yhi - ylo)/2 - buffer
    inner_radius = outer_radius - distance_between_membranes
    
    input_f.write('region membrane_inner cylinder x {} {} {} 0 {}\n'.format(
        yhi/2, zhi/2, inner_radius, xhi))
    input_f.write('region membrane_inner_tmp cylinder x {} {} {} 0 {} side out\n'.format(
        yhi/2, zhi/2, inner_radius, xhi))
    input_f.write('region membrane_outer cylinder x {} {} {} 0 {} side in\n\n'.format(
        yhi/2, zhi/2, outer_radius, xhi))
    
    input_f.write('region membrane intersect 2 membrane_inner_tmp membrane_outer\n\n')
    
    input_f.write('create_atoms 3 random {} {} membrane\n\n'.format(n_linkers_embedded, brn_seed))

    input_f.write('group chain1 type 1 2\n'.format())
    input_f.write('group chain_linkers type 2\n'.format())
    input_f.write('group embedded type 3\n\n'.format())

    input_f.write('pair_style lj/cut {}\n'.format(global_cutoff))
    input_f.write('pair_coeff * * {} {} {}\n'.format(0.0, 0.0, global_cutoff))
    input_f.write('pair_coeff 1 3 {} {} {}\n\n'.format(*memlinker_chain))

    for bst in bonds_styles:
        input_f.write('bond_style {}\n'.format(bst[0]))
        input_f.write('bond_coeff {} {} {}\n\n'.format(*bst[1:]))

    for ast in angle_styles:
        input_f.write('angle_style {}\n'.format(ast[0]))
        input_f.write('angle_coeff {} {} {}\n\n'.format(*ast[1:]))

    input_f.write('timestep {0:.10f}\n\n'.format(timestep))

    if (dump_minimization == 1):
        input_f.write('dump minimization all atom 1 dump.min.lammpstrj\n')
    
    input_f.write('minimize 0.0 1.0e-8 10000 10000\n\n')

    input_f.write('fix 1 all nve/limit 0.01\n\n')

    input_f.write(
    'fix wallchain1 chain1 wall/region membrane_inner lj93 {} {} {}\n'.format(*wall_chain))
    input_f.write('fix wallembed embedded wall/region membrane lj93 {} {} {}\n\n'.format(*wall_memlinker))
        
    if (linkers_on_filament == 1):
        input_f.write(
            'fix walllinkermem link_mem wall/region membrane lj93 {} {} {}\n\n'.format(*wall_linker))

    for i in range(n_linkers_filament):
        j = i + n_atoms + 1

        input_f.write('variable x{} equal x[{}]\n'.format(i+1, j))
        input_f.write('variable y{} equal y[{}]\n'.format(i+1, j))
        input_f.write('variable z{} equal z[{}]\n\n'.format(i+1, j))
    
    input_f.write('variable e1x equal x[1]\n')
    input_f.write('variable e1y equal y[1]\n')
    input_f.write('variable e1z equal z[1]\n\n')
    
    input_f.write('variable e2x equal x[{}]\n'.format(n_atoms))
    input_f.write('variable e2y equal y[{}]\n'.format(n_atoms))
    input_f.write('variable e2z equal z[{}]\n\n'.format(n_atoms))
    
    if (compute_CoM_MSD == 1):
    
        input_f.write('compute msdchain chain1 msd com yes average yes\n\n')
        
        input_f.write('variable comx equal c_msdchain[1]\n')
        input_f.write('variable comy equal c_msdchain[2]\n')
        input_f.write('variable comz equal c_msdchain[3]\n')
        input_f.write('variable comsq equal c_msdchain[4]\n\n')
    
    if (compute_CoM_pos == 1 or linkers_on_filament == 0):
        input_f.write('compute comchain chain1 com\n\n')
        
        input_f.write('variable comx equal c_comchain[1]\n')
        input_f.write('variable comy equal c_comchain[2]\n')
        input_f.write('variable comz equal c_comchain[3]\n\n')

    input_f.write('thermo_style custom step time temp etotal\n')
    input_f.write('thermo 10000\n\n')

    input_f.write('run {}\n\n'.format(10000))

    input_f.write('unfix 1\n')
    
    if (dump_minimization == 1):
        input_f.write('undump minimization\n\n')
    
    input_f.write('reset_timestep 0\n\n')
    
    input_f.write('variable tsteps equal time\n\n')

    if (linkers_on_filament == 1):
        input_f.write('fix link_pos all print {} "${{tsteps}} '.format(measure_distance_every))

        for i in range(n_linkers_filament):
            input_f.write('${{x{0}}} ${{y{0}}} ${{z{0}}} '.format(i+1))

        input_f.write('" file {} screen no\n\n'.format(link_pos_fname_str))
    
    input_f.write('fix e2e_pos all print {} "${{tsteps}} '.format(measure_distance_every))
    
    input_f.write('${e1x} ${e1y} ${e1z} ')
    input_f.write('${e2x} ${e2y} ${e2z} ')
    
    input_f.write('" file {} screen no\n\n'.format(end_pos_fname_str))
    
    if (compute_CoM_MSD == 1):      
        input_f.write('fix com_msd chain1 print {} "${{tsteps}} '.format(measure_distance_every))
        
        input_f.write('${comx} ${comy} ${comz} ${comsq}')
        
        input_f.write('" file {} screen no\n\n'.format(com_msd_fname_str))
    
    if (compute_CoM_pos == 1 or linkers_on_filament == 0):
        input_f.write('fix com_pos chain1 print {} "${{tsteps}} '.format(measure_distance_every))
        
        input_f.write('${comx} ${comy} ${comz}')
        
        input_f.write('" file {} screen no\n\n'.format(com_pos_fname_str))
    

    input_f.write('dump mydump all atom 1000 dump.lammpstrj\n\n')

    input_f.write(
        'fix 2 chain1 brownian {0} {1} gamma_t {2}\n'.format(*brownian))
    input_f.write('fix 3 embedded brownian {0} {1} gamma_t {2}\n\n'.format(*brownian_embedded))

    input_f.write('thermo 100000\n\n')

    input_f.write('run {}\n'.format(run_steps))
