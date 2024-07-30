# combine several train***.xyz
import os
import re
import shutil

lattice_count = 0
data_count = 0 

def merge_and_move_files(folder_path):
    # get all files
    files = os.listdir(folder_path)
    # path of train.xyz 
    output_file = os.path.join(folder_path, 'comb_input.data')

    # make train folder
    train_folder = os.path.join(folder_path, 'input_comb')
    os.makedirs(train_folder, exist_ok=True)

    with open(output_file, 'w') as merged_file:
        for file in files:
            if file.startswith('input') and file.endswith('.data') and file != 'input.data':
                # absolute path of file
                file_path = os.path.join(folder_path, file)
                # read file and write to train.
                with open(file_path, 'r') as f:
                    merged_file.write(f.read())


                # move file to train folder
                shutil.move(file_path, os.path.join(train_folder, file))


if __name__ == "__main__":
    folder_path = "."  # 同一フォルダ内で実行する場合
    merge_and_move_files(folder_path)

    #     # count dataset number
    # with open("comb_input.data", "r") as f:
    #     lattice_count += f.read().count("begin")
    #     lines = f.readlines()
    #     print(lines)
    #     for i in range(len(lines)):
    #         if 'comment' in lines[i]:
    #             data_count += 1
    #             print(data_count)
    #             lines[i] = "comment source_file_name=vasprun.xml structure_number={}".format(data_count)
    # print("The number of training datasets : " + str(lattice_count))
    # print("The combination and transfer of input.data files completed !")
    # os.rename("comb_input.data", "input.data")

    with open("comb_input.data", "r") as f:
        content = f.read()
        lattice_count += content.count("begin")
        lines = content.splitlines()

    # Print initial lines for debugging
    # print(lines)

    # Modify lines containing 'comment'
    for i in range(len(lines)):
        if 'comment' in lines[i]:
            data_count += 1
            lines[i] = "comment source_file_name=vasprun.xml structure_number={}".format(data_count)

    # Write the modified content back to the file
    with open("comb_input.data", "w") as f:
        f.write("\n".join(lines))

    # Output results
    print("The number of training datasets : " + str(lattice_count))
    print("The combination and transfer of input.data files completed !")

    # Rename the file
    os.rename("comb_input.data", "input.data")