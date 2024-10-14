# make IV mixed crystal with specific 6 membered rings !!
from ase import Atoms
from ase.build import bulk
from ase.io import write
from ase.neighborlist import NeighborList
import numpy as np
import random

# Dummy function to detect six-membered rings
# In practice, this should implement the real ring detection algorithm
num_rings = 0
mass = {
    "Si": 28.085,  # silicon
    "Ge": 72.63,   # Germanium
    "Sn": 118.71   # Tin
    } 

def replace_atom(atoms, index, symbol):
    atoms[index].symbol = symbol

def get_input_element(prompt):
    element = input(prompt).strip()
    if element not in ['Si', 'Ge', 'Sn']:
        raise ValueError("Please select one of Si, Ge, or Sn.")
    return element

def get_input_cell_size():
    cell_size = input("Enter the cell size as three integers (e.g., '4 4 10'): ")
    return tuple(map(int, cell_size.split()))

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
    #print(path)
    if depth == max_depth:
        # Check if the last atom is connected to the first atom to close the ring
        neighbors = nl.get_neighbors(atom_index)[0]
        if path[0] in neighbors:
            visited.add(atom_index)
            path.append(atom_index)
            # Check if the ring is homogeneous (Si or Ge only)
            if is_homogeneous_ring(path, structure):
                ring_type = classify_ring(path, structure)
                return [(path[:], ring_type)]  # Return the ring and its type
        return []  # No valid ring found

    # Add the current atom to the visited set and path
    visited.add(atom_index)
    path.append(atom_index)
    
    # Get neighbors of the current atom
    neighbors = nl.get_neighbors(atom_index)[0]
    rings = []
    
    # Visit each neighbor if it's not already in the visited set and not creating loops prematurely
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
    
    for atom_index in range(len(structure)):
        visited = set()  # To track visited atoms for each DFS search
        path = []  # The current path in DFS
        detected_rings = dfs(atom_index, structure, nl, visited, path, 1)  # Start DFS from this atom
        
        for ring in detected_rings[0]:
            if ring != "Si" and ring != "Ge" and ring != "Sn":
                #print(ring)
                rings.append((list(ring)))

    
    # Optionally, remove duplicate rings (since rings can be found starting from different atoms)
    #print(rings)
    return rings

def find_center_atom(atoms):
    center = np.mean(atoms.get_positions(), axis=0)
    distances = np.linalg.norm(atoms.get_positions() - center, axis=1)
    print(f"The ID of the center of crystal : {np.argmin(distances)}")
    return np.argmin(distances)

def main():
    # ① Prompt the user to select an atom type from Si, Ge, Sn
    main_element = get_input_element("Select the main crystal element (Si, Ge, Sn): ")
    secondary_elements = [el for el in ['Si', 'Ge', 'Sn'] if el != main_element]
    ring_element = get_input_element(f"Select the element for six-membered rings from {secondary_elements}: ")

    # ② Create a diamond structure and prompt the user for cell expansion size
    cell_size = get_input_cell_size()
    
    # Build diamond structure (Si, Ge, Sn are supported as examples)
    crystal = bulk(main_element, 'diamond', a=5.431, cubic=True)  # a=5.431 is the lattice constant for Si
    crystal = crystal * cell_size  # Expand the cell

    # Prompt the user for the required number of six-membered rings
    required_rings = int(input("Enter the required number of six-membered rings: "))

    # ③ Get the central atom's ID and detect all six-membered rings
    center_atom_id = find_center_atom(crystal) # Assume the central atom is the midpoint of the entire cell
    
    # six-membered analysis and output
    radii = 2.4  # Approximate cutoff for Si-Ge bonds
    nl = NeighborList([radii / 2] * len(crystal), skin=0.3, bothways=True, self_interaction=False)
    nl.update(crystal)
    all_rings = find_six_membered_rings(crystal, nl)
    #print(f"Number of {main_element} six-membered rings: {len(all_rings)}")

    # Increase rings candidate
    substitution_candidates = []
    candidate_rings = [ring for ring in all_rings if center_atom_id in ring]

    select_rings = []
    if len(candidate_rings) >= required_rings:
        select_rings = random.sample(candidate_rings, required_rings)
        for ring in select_rings:
            #print(ring)
            substitution_candidates.extend(ring)
        substitution_candidates = list(set(substitution_candidates))
        #print(substitution_candidates)

    while len(candidate_rings) < required_rings:
        for ring in candidate_rings:
            substitution_candidates.extend(ring)
        substitution_candidates = list(set(substitution_candidates))
        for idx in substitution_candidates:
            candidate_rings.extend([ring for ring in all_rings if int(idx) in ring])
            candidate_rings = list(set(tuple(ring) for ring in candidate_rings))
            if len(candidate_rings) > required_rings:
                #print(len(candidate_rings))
                select_rings = candidate_rings[:required_rings]
                for ring in select_rings:
                    substitution_candidates.extend(ring)
                substitution_candidates = list(set(substitution_candidates))
                break

    
    # replace atom
    for idx in substitution_candidates:
            replace_atom(crystal, idx, str(ring_element))

    print(f"All atom : {len(crystal)}")
    print(f"{len(substitution_candidates)} {main_element} atoms -> {ring_element} atoms transformed.")

    # ⑤ Find rings that share two or more atoms with the substitution candidate list and add them
    #for _ in range(required_rings - 1):
        #new_candidates = []
        #for ring in all_rings:
            #if len(set(ring).intersection(substitution_candidates)) >= 2:
                #new_candidates.extend(ring)
        #substitution_candidates.extend(new_candidates)
        #substitution_candidates = list(set(substitution_candidates))  # Remove duplicates

    #print(substitution_candidates)

    # replace atoms
    #for idx in substitution_candidates:
        #replace_atom(crystal, idx, str(ring_element))

    # ⑥ Write the structure in LAMMPS data format
    # Add Masses block to lmp file
    # Reading the file contents
    crystal.write('cycro-cluster.lmp', format='lammps-data')
    file_path = 'cycro-cluster.lmp'
    with open(file_path, 'r') as file:
        file_content = file.readlines()

    # Lines to be inserted after line 8
    lines_to_insert = [
    '\nMasses\n',
    '\n',
    f'1 {mass[ring_element]}\n',
    f'2 {mass[main_element]}\n',
    '\n'
    ]

    # Inserting the lines at the 9th position (after index 8)
    updated_content = file_content[:8] + lines_to_insert + file_content[8:]

    # Saving the modified file
    with open(file_path, 'w') as modified_file:
        modified_file.writelines(updated_content)

        print(f'Structure saved in {file_path}')

if __name__ == '__main__':
    main()
