import os

# Function to convert input.data to train.xyz
def convert_input_data_to_train_xyz(input_data_path, output_train_xyz_path):
    with open(input_data_path, 'r') as input_file:
        lines = input_file.readlines()
    
    with open(output_train_xyz_path, 'w') as output_file:
        block = []
        lattice = []
        atoms = []
        energies = []
        comment = ""
        
        for line in lines:
            if line.startswith('begin'):
                block = []
                lattice = []
                atoms = []
                energies = []
            elif line.startswith('lattice'):
                lattice.append(line.split()[1:])
            elif line.startswith('atom'):
                atoms.append(line.split()[1:])
            elif line.startswith('energy'):
                energies.append(line.split()[1])
            elif line.startswith('comment'):
                comment = line.strip().split(maxsplit=1)[1]
            elif line.startswith('end'):
                # Write to output file in train.xyz format
                num_atoms = len(atoms)
                output_file.write(f'{num_atoms}\n')
                
                # Construct lattice string
                lattice_str = " ".join(sum(lattice, []))
                
                # Construct energy strings
                energy_str = energies[0]
                free_energy_str = energies[0]
                
                output_file.write(f'Lattice="{lattice_str}" Properties=species:S:1:pos:R:3:forces:R:3 energy={energy_str} free_energy={free_energy_str} pbc="T T T"\n')
                
                for atom in atoms:
                    output_file.write(f'{atom[3]} {atom[0]} {atom[1]} {atom[2]} {atom[5]} {atom[6]} {atom[7]}\n')

# Function to convert train.xyz to input.data
def convert_train_xyz_to_input_data(train_xyz_path, output_input_data_path):
    with open(train_xyz_path, 'r') as input_file:
        lines = input_file.readlines()
    
    with open(output_input_data_path, 'w') as output_file:
        i = 0
        while i < len(lines):
            if lines[i].strip().isdigit():
                num_atoms = int(lines[i].strip())
                i += 1
                properties_line = lines[i].strip()
                lattice_str = properties_line.split('Lattice="')[1].split('"')[0]
                lattice = lattice_str.split()
                
                energy_str = properties_line.split('energy=')[1].split()[0]
                free_energy_str = properties_line.split('free_energy=')[1].split()[0]
                
                output_file.write('begin\n')
                output_file.write(f'lattice {lattice[0]} {lattice[1]} {lattice[2]}\n')
                output_file.write(f'lattice {lattice[3]} {lattice[4]} {lattice[5]}\n')
                output_file.write(f'lattice {lattice[6]} {lattice[7]} {lattice[8]}\n')
                
                for j in range(num_atoms):
                    i += 1
                    atom_line = lines[i].strip().split()
                    species = atom_line[0]
                    pos_x, pos_y, pos_z = atom_line[1], atom_line[2], atom_line[3]
                    force_x, force_y, force_z = atom_line[4], atom_line[5], atom_line[6]
                    output_file.write(f'atom {pos_x} {pos_y} {pos_z} {species} 0.0 0.0 {force_x} {force_y} {force_z}\n')
                
                output_file.write(f'energy {energy_str}\n')
                output_file.write('comment Converted from train.xyz\n')
                output_file.write('end\n')
            
            i += 1

# Main program
def main():
    # Get the list of files in the current directory
    files = [f for f in os.listdir('.') if os.path.isfile(f) and (f.endswith('.xyz') or f.endswith('.data'))]
    
    if not files:
        print("No files to convert found.")
        return
    
    # Display the list of files
    print("Select a file to convert:")
    for i, file in enumerate(files):
        print(f"{i + 1}: {file}")
    
    # Select a file
    while True:
        try:
            choice = int(input("Enter the file number: ")) - 1
            if 0 <= choice < len(files):
                selected_file = files[choice]
                break
            else:
                print("Invalid number. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")
    
    # Perform conversion
    if selected_file.endswith('.data'):
        output_file = selected_file.replace('.data', '_converted.xyz')
        convert_input_data_to_train_xyz(selected_file, output_file)
    elif selected_file.endswith('.xyz'):
        output_file = selected_file.replace('.xyz', '_converted.data')
        convert_train_xyz_to_input_data(selected_file, output_file)
    
    print(f"Conversion complete: {output_file}")

if __name__ == "__main__":
    main()
