from ase.io import read
from ase.build import bulk
from ase.visualize import view
from ase.neighborlist import NeighborList
import numpy as np
import random

# neighbor listで結合管理？

def get_input_cell_size():
    # Ask user for the cell size
    cell_size = input("Enter the cell size as three integers (e.g., '4 4 10'): ")
    return tuple(map(int, cell_size.split()))

def replace_randomly_with_ge(atoms, num_ge, seed):
    random.seed(seed)
    indices = list(range(len(atoms)))
    ge_indices = random.sample(indices, num_ge)
    for idx in ge_indices:
        atoms[idx].symbol = 'Ge'
    return atoms

def create_structure():
    # Step 1: Create the initial Si diamond structure
    cell_size = get_input_cell_size()
    si_structure = bulk('Si', 'diamond', a=5.43, cubic=True) * cell_size
    
    num_atoms = len(si_structure)
    num_ge = num_atoms // 2  # Half of the atoms should be Ge
    # Replace atoms randomly, with seed input
    seed = int(input("Enter a seed value for random replacement: "))
    si_ge_structure = replace_randomly_with_ge(si_structure, num_ge, seed)

    return si_ge_structure

def is_homogeneous_ring(ring, structure):
    """
    Check if all atoms in the ring are either Si or Ge.
    
    Parameters:
    - ring: A list of atom indices representing the ring.
    - structure: The atomic structure object from ASE.
    
    Returns:
    - True if the ring is composed of only Si or only Ge atoms.
    """
    symbols = [structure[i].symbol for i in ring]
    return all(s == 'Si' for s in symbols) or all(s == 'Ge' for s in symbols)

def classify_ring(ring, structure):
    """
    Classify whether a ring is made up of Si or Ge atoms.
    
    Parameters:
    - ring: A list of atom indices representing the ring.
    - structure: The atomic structure object from ASE.
    
    Returns:
    - A string indicating if the ring is 'Si' or 'Ge', or None if it's mixed (shouldn't happen).
    """
    symbols = [structure[i].symbol for i in ring]
    if all(s == 'Si' for s in symbols):
        return 'Si'
    elif all(s == 'Ge' for s in symbols):
        return 'Ge'
    else:
        return None  # Shouldn't happen, as mixed rings are filtered out

def dfs(atom_index, structure, nl, visited, path, depth, max_depth=6):
    """
    Perform depth-first search (DFS) to detect six-membered rings, filtering out mixed Si/Ge rings.
    
    Parameters:
    - atom_index: The current atom in the DFS path.
    - structure: The atomic structure object from ASE.
    - nl: NeighborList for finding neighbors of each atom.
    - visited: A set of atoms already visited in this path.
    - path: The current path of atoms being explored.
    - depth: Current depth in the DFS (number of atoms in the path).
    - max_depth: The depth at which we expect a ring (6 for six-membered rings).
    
    Returns:
    - A list of tuples where each tuple contains the atom indices of a ring and its type ('Si' or 'Ge').
    """
    if depth > max_depth:
        return []
    
    # Add the current atom to the visited set and path
    visited.add(atom_index)
    path.append(atom_index)
    
    # If we've reached a depth of 6, check if we have a valid ring
    if depth == max_depth:
        # Check if the last atom in the path is connected to the first
        neighbors = nl.get_neighbors(atom_index)[0]
        if path[0] in neighbors:
            # Check if the ring is homogeneous (Si or Ge only)
            if is_homogeneous_ring(path, structure):
                ring_type = classify_ring(path, structure)
                return [(path[:], ring_type)]  # Return the ring and its type
    
    # Get neighbors of the current atom
    neighbors = nl.get_neighbors(atom_index)[0]
    rings = []
    
    # Visit each neighbor if it's not already in the visited set
    for neighbor in neighbors:
        if neighbor not in visited:
            rings.extend(dfs(neighbor, structure, nl, visited, path, depth + 1, max_depth))
    
    # Backtrack: remove the current atom from the path and visited set
    path.pop()
    visited.remove(atom_index)
    
    return rings

def find_six_membered_rings(structure, nl):
    """
    Find six-membered homogeneous rings (Si or Ge only) in the structure using DFS.
    
    Parameters:
    - structure: The atomic structure object from ASE.
    - nl: NeighborList for finding neighbors.
    
    Returns:
    - A list of tuples where each tuple contains a six-membered homogeneous ring's atom indices and its type ('Si' or 'Ge').
    """
    rings = []
    si_count = 0
    ge_count = 0
    
    for atom_index in range(len(structure)):
        visited = set()  # To track visited atoms for each DFS search
        path = []  # The current path in DFS
        detected_rings = dfs(atom_index, structure, nl, visited, path, 1)  # Start DFS from this atom
        
        # Count Si and Ge rings
        for ring, ring_type in detected_rings:
            rings.append((ring, ring_type))
            if ring_type == 'Si':
                si_count += 1
            elif ring_type == 'Ge':
                ge_count += 1
    
    # Optionally, remove duplicate rings (since rings can be found starting from different atoms)
    unique_rings = [list(t) for t in {tuple(sorted(ring[0])): ring[1] for ring in rings}.items()]
    
    return unique_rings, si_count, ge_count


# Generate the structure and view it
structure = create_structure()
view(structure)  # Optional, to visualize the structure
structure.write('cycro-SiGe.lmp', format='lammps-data')

# Load the provided file and inspect the first few lines to identify the right position for the insertion.
file_path = 'SiGe-cycro.lmp'

# six-membered analysis and output
radii = 2.4  # Approximate cutoff for Si-Ge bonds
nl = NeighborList([radii / 2] * len(structure), skin=0.3, bothways=True, self_interaction=False)
nl.update(structure)
#rings = find_six_membered_rings(structure, nl)
rings, si_count, ge_count = find_six_membered_rings(structure, nl)

for ring, ring_type in rings:
    print(f"Ring of type {ring_type} with atoms: {ring}")
print(f"Number of homogeneous six-membered rings found: {len(rings)}")
si_ring = 0
ge_ring = 0
for i in range(len(rings)):
    if rings[i][1] == 'Si':
        si_ring += 1
    elif rings[i][1] == 'Ge':
        ge_ring += 1
print(f"Number of Si six-membered rings: {si_ring}")
print(f"Number of Ge six-membered rings: {ge_ring}")

# Add Masses block to lmp file
# Reading the file contents
with open(file_path, 'r') as file:
    file_content = file.readlines()

# Lines to be inserted after line 8
lines_to_insert = [
    'Masses\n',
    '\n',
    '1 28.0855\n',
    '2 58.6934\n',
    '\n'
]

# Inserting the lines at the 9th position (after index 8)
updated_content = file_content[:8] + lines_to_insert + file_content[8:]

# Saving the modified file
with open(file_path, 'w') as modified_file:
    modified_file.writelines(updated_content)

