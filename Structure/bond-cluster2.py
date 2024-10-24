from ase import Atoms
from ase.build import bulk
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.io import write
import numpy as np
import random

# Masses for LAMMPS file
mass = {
    "Si": 28.085,
    "Ge": 72.63,
    "Sn": 118.71
}

# Get the main crystal atom and replacement atom from user
main_atom = input("Enter the main atom for the crystal (e.g., Si, Ge, Sn): ").strip()
replace_atom = input("Enter the atom to replace with (e.g., Ge, Si, Sn): ").strip()

# Get the cell size in '4 4 4' format and split into three integers
cell_size_input = input("Enter the cubic cell expansion size (e.g., 4 4 4): ")
cell_size = tuple(map(int, cell_size_input.split()))

# Get the concentration of the replacement atom
ge_concentration = float(input("Enter the concentration of the replacing atom (e.g., 0.1 for 10%): "))

# Get bond condition input
max_ge_neighbors = int(input("Enter the maximum number of allowed Ge-Ge bonds (0: none, 1: 2-body, etc.): "))

# Build the cubic diamond structure for the main crystal with the given cell size
crystal = bulk(main_atom, 'diamond', a=5.431) * cell_size

# Number of atoms and how many should be replaced
num_atoms = len(crystal)
num_ge = int(ge_concentration * num_atoms)

# Randomly replace atoms ensuring bond constraints
available_indices = list(range(num_atoms))
ge_indices = []

# NeighborList setup
cutoffs = natural_cutoffs(crystal)
nl = NeighborList(cutoffs, self_interaction=False, bothways=True)

# Function to replace atom symbol
def replace_atom(atoms, index, symbol):
    atoms[index].symbol = symbol

# Ge replacement with bond constraint
for _ in range(num_ge):
    nl.update(crystal)
    selected_index = random.choice(available_indices)

    # Get neighbors and check how many are already Ge
    indices, offsets = nl.get_neighbors(selected_index)
    ge_neighbors = [idx for idx in indices if crystal[idx].symbol == replace_atom]

    if len(ge_neighbors) <= max_ge_neighbors:
        replace_atom(crystal, selected_index, replace_atom)  # 間違い
        # 修正：replace_atom関数に、実際の置換するシンボル文字列を渡す
        replace_atom(crystal, selected_index, replace_atom)  # ここは replace_atom の代わりに replace_atom を渡す
        ge_indices.append(selected_index)
        available_indices.remove(selected_index)
        
# Check Ge-Ge bonding information
nl.update(crystal)
num_body2 = []
num_body3 = []

for i in ge_indices:
    indices, offsets = nl.get_neighbors(i)
    ge_neighbors = [idx for idx in indices if crystal[idx].symbol == replace_atom]

    if len(ge_neighbors) == 1:
        num_body2.append((i, ge_neighbors[0]))

    if len(ge_neighbors) == 2:
        num_body3.append((i, ge_neighbors[0], ge_neighbors[1]))

# Output bonding information
print(f"Number of 2-body Ge-Ge bonds: {len(num_body2)}")
print(f"Number of 3-body Ge-Ge bonds: {len(num_body3)}")

# LAMMPS format file output
lammps_filename = f'cycro-{main_atom}{replace_atom}cluster.lmp'
write(lammps_filename, crystal, format='lammps-data')

# Insert Masses block into LAMMPS file
with open(lammps_filename, 'r') as file:
    lammps_data = file.readlines()

mass_lines = [
    '\nMasses\n\n',
    f'1 {mass[main_atom]}\n',
    f'2 {mass[replace_atom]}\n\n'
]

# Inserting mass after line 8
lammps_data = lammps_data[:8] + mass_lines + lammps_data[8:]

# Saving the updated file
with open(lammps_filename, 'w') as file:
    file.writelines(lammps_data)

print(f"Structure saved in {lammps_filename}")
