from genericpath import isfile
import os
import re
import shutil
import argparse


# change in case of PC or number of parallel
path_lmp = "lmp_allegro"

# get in. file
input_file_list = []
files_and_dirs = os.listdir('.')
# serach in. file
for item in files_and_dirs:
    if os.path.isfile(item) and item.startswith('in.'):
        input_file_list.append(item)
        print(item)
print("Select input file by enter number")
for i in range(len(input_file_list)):
    print(i + " : " + input_file_list[i])
index_input = int(input("Enter the number : "))
input_file = input_file_list[index_input]

# get .lmp file
for item in files_and_dirs:
    if os.path.isfile(item) and item.endswith('.lmp'):
        lmp_file = item

# get .lmp file
for item in files_and_dirs:
    if os.path.isfile(item) and item.endswith('.pth'):
        pth_file = item

# scaling factor for 2nd line in POSCAR
start = 0.9
end = 1.1
step = 0.05 
scaler = [float(x/100) for x in range(int(start * 100), int(end * 100) + 1, int(step * 100))]
dir_num = 0


parser = argparse.ArgumentParser(description="input two argument")
parser.add_argument("-d", required=True, help="Please input directory name")
args = parser.parse_args()

path = os.getcwd()
dirlist = os.listdir(path)
dir_list = []
for i in dirlist:
    if os.path.isdir(i):
        dir_list.append(i)
while True:
    flag = 0
    for i in dir_list:
        if re.search(args.d, i):
            print("The directory named {} already exists".format(args.d))
            flag = 1
            break
    if flag == 1:
        exit()
    else:
        break


for calc_num in scaler:
    dir_num+=1
    dir_name = args.d + "_" +"{:.02f}".format(calc_num)+"_"+str(dir_num)
    os.makedirs(dir_name)
    
    if dir_num == 1:
        # copy files
        shutil.copy(input_file, dir_name)
        shutil.copy(lmp_file, dir_name)
        if os.path.isfile(pth_file):    
            shutil.copy(pth_file, dir_name)        
        
        # Define the new line to be added
        new_line = "change_box all x scale "+ str(calc_num) +" y scale "+ str(calc_num) + "z scale "+ str(calc_num) +"\n"

        # Read the file
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # List to store the new content
        new_lines = []

        # Process the file content and insert the new line
        for i in range(len(lines)):
            if lines[i].startswith('run'):
                # Add the new line above the line starting with "run"
                new_lines.append(new_line)
            # Add the current line to the new list
            new_lines.append(lines[i])

        # Write the new content back to the file
        with open(input_file, 'w') as f:
            f.writelines(new_lines)

        # execute
        os.system(path_lmp + " -in " + input_file)
        os.chdir("..")
        
    else:
        # dir_name_b = args.d + "_" +  str(float((int(100*calc_num)-int(100*step))/100)) +"_"+str(dir_num-1)
        dir_name_b = args.d + "_" +  "{:.02f}".format(scaler[dir_num-2]) +"_"+str(dir_num-1)
        
        # copy files
        shutil.copy(input_file, dir_name)
        shutil.copy(lmp_file, dir_name)
        if os.path.isfile(pth_file):    
            shutil.copy(pth_file, dir_name) 

        
        # write scaling factor
        # Define the new line to be added
        new_line = "change_box all x scale "+ str(calc_num) +" y scale "+ str(calc_num) + "z scale "+ str(calc_num) +"\n"

        # Read the file
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # List to store the new content
        new_lines = []

        # Process the file content and insert the new line
        for i in range(len(lines)):
            if lines[i].startswith('run'):
                # Add the new line above the line starting with "run"
                new_lines.append(new_line)
            # Add the current line to the new list
            new_lines.append(lines[i])

        # Write the new content back to the file
        with open(input_file, 'w') as f:
            f.writelines(new_lines)

        # execute
        os.system(path_lmp + " -in " + input_file)
        os.chdir("..")

