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
                    output_file.write(f'{atom[3]} {atom[0]} {atom[1]} {atom[2]} {atom[6]} {atom[7]} {atom[8]}\n')
            
convert_input_data_to_train_xyz('input100_11.data', 'train-11.xyz')
