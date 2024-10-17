import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.neighborlist import neighbor_list
from ase.io import write
import random

# Atomic masses dictionary
mass = {
    "Si": 28.085,  # silicon
    "Ge": 72.63,   # germanium
    "Sn": 118.71   # tin
}

# Replace a specific atom in the atoms object
def replace_atom(atoms, index, symbol):
    atoms[index].symbol = symbol

# Get input for the element (Si, Ge, Sn)
def get_input_element(prompt):
    element = input(prompt).strip()
    if element not in ['Si', 'Ge', 'Sn']:
        raise ValueError("Please select one of Si, Ge, or Sn.")
    return element

# Get input for the cell size
def get_input_cell_size():
    cell_size = input("Enter the cell size as three integers (e.g., '4 4 10'): ")
    return tuple(map(int, cell_size.split()))

# Get replacement concentration and bonding limit (2-body or 3-body)
def get_replacement_settings():
    concentration = float(input("Enter the replacement concentration (e.g., 0.2 for 20%): "))
    bond_limit = int(input("Enter the bonding limit (2 for two-body or 3 for three-body): "))
    if bond_limit not in [2, 3]:
        raise ValueError("Please select either 2 or 3 for bonding limit.")
    max_iterations = int(input("Enter the maximum number of iterations: "))
    return concentration, bond_limit, max_iterations

# Find the center atom in the structure
def find_center_atom(atoms):
    center = np.mean(atoms.get_positions(), axis=0)
    distances = np.linalg.norm(atoms.get_positions() - center, axis=1)
    print(f"The ID of the center of crystal : {np.argmin(distances)}")
    return np.argmin(distances)

# Create the crystal structure based on user input
def create_crystal():
    base_element = get_input_element("Enter the base element (Si, Ge, Sn): ")
    replace_element = get_input_element("Enter the replacement element (Si, Ge, Sn, but different from base): ")
    if base_element == replace_element:
        raise ValueError("Replacement element must be different from base element.")
    
    cell_size = get_input_cell_size()
    lattice_constant = 5.43  # Default lattice constant for Si/Ge/Sn (can be adjusted)
    atoms = bulk(base_element, 'diamond', a=lattice_constant, cubic=True).repeat(cell_size)
    
    return atoms, base_element, replace_element

# Perform atom replacement based on the given concentration
def perform_replacement(atoms, base_element, replace_element, concentration, avoid_indices):
    num_atoms = len(atoms)
    num_replace_atoms = int(concentration * num_atoms)  # Replace based on concentration
    available_indices = [i for i in range(num_atoms) if i not in avoid_indices]
    
    # Randomly select atoms to replace, avoiding those in avoid_indices
    replace_indices = random.sample(available_indices, num_replace_atoms)
    
    for idx in replace_indices:
        replace_atom(atoms, idx, replace_element)
    
    return replace_indices

# Detect two-body and three-body bonds only between replaced atoms
def detect_bonds(atoms, replace_indices, bond_limit):
    cutoff = 2.8  # Cutoff distance for bonds (Å)
    i, j, dists = neighbor_list('ijd', atoms, cutoff=cutoff)
    
    # Two-body bond detection (only between replaced atoms)
    two_body_bonds = []
    three_body_bonds = []
    
    for idx1, idx2 in zip(i, j):
        if idx1 in replace_indices and idx2 in replace_indices:
            if idx1 < idx2:  # Avoid duplicates
                two_body_bonds.append((idx1, idx2))
    
    # Three-body bond detection (if allowed by bond_limit)
    if bond_limit == 3:
        for idx1 in replace_indices:
            neighbors = [idx2 for idx2 in j[i == idx1] if idx2 in replace_indices]
            if len(neighbors) >= 2:
                for idx2 in neighbors:
                    for idx3 in neighbors:
                        if idx2 < idx3:
                            three_body_bonds.append((idx1, idx2, idx3))
    
    print(f"Two-body bonds between replaced atoms: {two_body_bonds}")
    print(f"Three-body bonds between replaced atoms: {three_body_bonds}")
    
    return two_body_bonds, three_body_bonds

# Get the indices to avoid in the next replacement step
def get_avoid_indices(two_body_bonds, three_body_bonds, bond_limit):
    avoid_indices = set()
    for bond in two_body_bonds:
        avoid_indices.update(bond)
    if bond_limit == 3:
        for bond in three_body_bonds:
            avoid_indices.update(bond)
    return list(avoid_indices)

# Replace atoms until no bonds are left, or max iterations reached
def refine_structure(atoms, base_element, replace_element, concentration, bond_limit, max_iterations):
    replace_indices = perform_replacement(atoms, base_element, replace_element, concentration, [])
    
    for iteration in range(max_iterations):
        two_body_bonds, three_body_bonds = detect_bonds(atoms, replace_indices, bond_limit)
        
        if bond_limit == 2 and len(two_body_bonds) == 0 and len(three_body_bonds) == 0:
            print(f"No more two-body or three-body bonds left after {iteration} iterations.")
            break
        elif bond_limit == 3 and len(three_body_bonds) == 0:
            print(f"No more three-body bonds left after {iteration} iterations.")
            break
        
        # Get the indices to avoid in the next replacement step
        avoid_indices = get_avoid_indices(two_body_bonds, three_body_bonds, bond_limit)
        
        # Remove existing bonds by replacing atoms back to base and avoid certain atoms
        if bond_limit == 2:
            for idx1, idx2 in two_body_bonds:
                replace_atom(atoms, idx2, base_element)  # Replace bonded atom back to base
            replace_indices = perform_replacement(atoms, base_element, replace_element, concentration, avoid_indices)
        
        elif bond_limit == 3:
            for idx1, idx2, idx3 in three_body_bonds:
                replace_atom(atoms, idx2, base_element)  # Replace bonded atom back to base
                replace_atom(atoms, idx3, base_element)
            replace_indices = perform_replacement(atoms, base_element, replace_element, concentration, avoid_indices)
    
    if iteration == max_iterations - 1:
        print(f"Reached maximum iterations ({max_iterations}) without fully removing bonds.")

# Output the structure to LAMMPS data format with mass info
def output_structure(atoms, base_element, replace_element, bond_limit):
    with open(f'{bond_limit}bond-{base_element}{replace_element}.lmp', 'w') as f:
        write(f, atoms, format='lammps-data')
        f.write("\nMasses\n\n")
        f.write(f"1 {mass[base_element]}\n")
        f.write(f"2 {mass[replace_element]}\n")
    print(f"The structure with mass info has been saved as {bond_limit}bond-{base_element}{replace_element}.lmp.")

# Main execution
def main():
    atoms, base_element, replace_element = create_crystal()
    concentration, bond_limit, max_iterations = get_replacement_settings()
    
    # Refine structure to remove unwanted bonds
    refine_structure(atoms, base_element, replace_element, concentration, bond_limit, max_iterations)
    
    # Output the final structure
    output_structure(atoms, base_element, replace_element, bond_limit)

main()
