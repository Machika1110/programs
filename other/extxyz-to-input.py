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
                output_file.write('comment Converted from train.xyz\n')
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
                output_file.write('charge 0.0\n')
                output_file.write('end\n')
            
            i += 1
            
convert_train_xyz_to_input_data('train100.xyz', 'input-100.data')
