# combine several train***.xyz
import os
import shutil

def merge_and_move_files(folder_path):
    # get all files
    files = os.listdir(folder_path)
    # path of train.xyz 
    output_file = os.path.join(folder_path, 'comb_train.xyz')

    # maake train folder
    train_folder = os.path.join(folder_path, 'train')
    os.makedirs(train_folder, exist_ok=True)

    lattice_count = 0 

    with open(output_file, 'w') as merged_file:
        for file in files:
            if file.startswith('train') and file.endswith('.xyz') and file != 'comb_train.xyz':
                # absolute path of file
                file_path = os.path.join(folder_path, file)
                # read file and write to train.
                with open(file_path, 'r') as f:
                    merged_file.write(f.read())


                # move file to train folder
                shutil.move(file_path, os.path.join(train_folder, file))
    # count dataset number
    with open("comb_train.xyz", "r") as f:
        lattice_count += f.read().count("Lattice")
    print("The number of training datasets : " + str(lattice_count))
    print("The combination and transfer of train files completed !")
    os.rename("comb_train.xyz", "train.xyz")


if __name__ == "__main__":
    folder_path = "."  # 同一フォルダ内で実行する場合
    merge_and_move_files(folder_path)