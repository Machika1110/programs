from ase.io import read, write
from ase.build import bulk
from ase.visualize import view
from ase.neighborlist import NeighborList
import numpy as np
import random
from collections import deque
import os

# Helper functions
def get_input_element(prompt):
    element = input(prompt).strip()
    if element not in ['Si', 'Ge', 'Sn']:
        raise ValueError("Please select one of Si, Ge, or Sn.")
    return element

def get_input_cell_size():
    cell_size = input("Enter the cell size as three integers (e.g., '4 4 10'): ")
    return tuple(map(int, cell_size.split()))

def get_ring_count():
    return int(input("How many six-membered rings do you want to form? "))

# Replace atoms in bulk structure
def replace_atom(atoms, index, symbol):
    atoms[index].symbol = symbol

# Create crystal structure with given dimensions and element
def create_structure(element, cell_size):
    structure = bulk(element, 'diamond', a=5.43, cubic=True) * cell_size
    structure.write('seed-crystal.lmp', format='lammps-data')
    return structure

# Find the atom closest to the center of the structure
def find_center_atom(atoms):
    center = np.mean(atoms.get_positions(), axis=0)
    distances = np.linalg.norm(atoms.get_positions() - center, axis=1)
    return np.argmin(distances)

# Depth-First Search to form rings
def dfs_replace_to_form_ring(start_index, atoms, neighbor_list, target_element, formed_rings, target_ring_size=6):
    """
    This function uses a depth-first search (DFS) to explore paths from a given atom, aiming to form a 
    ring of a specified size (e.g., 6-membered ring). When such a ring is found, the atoms involved in 
    the ring are replaced by the target element.
    
    Parameters:
    - start_index: Index of the starting atom.
    - atoms: The atomic structure from ASE.
    - neighbor_list: NeighborList from ASE to track bonded neighbors.
    - target_element: The element to replace the atoms in the ring.
    - target_ring_size: Size of the ring to form (default is 6).
    
    Returns:
    - True if a ring was formed and atoms were replaced, False otherwise.
    """
    visited = set()
    stack = [(start_index, [])]  # (current atom index, current path)

    # Perform DFS to find a ring
    while stack:
        current, path = stack.pop()
        # Skip already visited atoms
        if current in visited:
            continue
        visited.add(current)
        path.append(current)
        #print(len(path))
        # If we reach a path of length target_ring_size, check if it's a ring
        if len(path) == target_ring_size * (formed_rings+1):
            # Check if the current atom is a neighbor of the start atom, forming a ring
            neighbors, _ = neighbor_list.get_neighbors(current)
            #if start_index in neighbors:
            if True:
                # We have found a ring; replace atoms in this path
                for idx in path:
                    replace_atom(atoms, idx, target_element)
                return True

        # Explore the neighbors of the current atom
        neighbors, _ = neighbor_list.get_neighbors(current)
        for neighbor in neighbors:
            if neighbor not in visited:
                stack.append((neighbor, path.copy()))  # Continue DFS with this neighbor

    return False

# Main process
def main():
    # Step 1: Get inputs from user
    main_element = get_input_element("Select the main crystal element (Si, Ge, Sn): ")
    secondary_elements = [el for el in ['Si', 'Ge', 'Sn'] if el != main_element]
    ring_element = get_input_element(f"Select the element for six-membered rings from {secondary_elements}: ")
    
    cell_size = get_input_cell_size()
    num_rings = get_ring_count()

    # Step 2: Create the main structure
    atoms = create_structure(main_element, cell_size)

    # Step 3: Find the center atom and replace with ring element
    center_atom_idx = find_center_atom(atoms)
    #print(center_atom_idx)
    replace_atom(atoms, center_atom_idx, ring_element)

    # Step 4: Setup neighbor list for DFS
    cutoffs = [2.4] * len(atoms)  # Approximate cutoff for covalent bonding in Si, Ge, Sn
    neighbor_list = NeighborList(cutoffs, skin=0, sorted=False, self_interaction=False)
    #print(neighbor_list)
    neighbor_list.update(atoms)

    # Step 5: Form rings iteratively
    formed_rings = 0
    os.mkdir("traj")
    while formed_rings < num_rings:
        success = dfs_replace_to_form_ring(center_atom_idx, atoms, neighbor_list, ring_element, formed_rings)
        atoms.write(f'traj/cycro-cluster-{formed_rings}.lmp', format='lammps-data')
        if not success:
            print(f"Failed to form more six-membered rings after {formed_rings} rings.")
            break
        formed_rings += 1

        # Save trajectory after each ring is formed
        # write(f"ring_{formed_rings}.traj", atoms)
        print(f"Formed {formed_rings} ring(s). Trajectory saved as 'ring_{formed_rings}.traj'.")
    
    atoms.write('cycro-cluster.lmp', format='lammps-data')
    print("Process completed.")


# Entry point
if __name__ == '__main__':
    main()

