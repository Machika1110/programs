# this program can be used for lectangier (α=β=γ=90°)
# "l_range" decides range of not moving
while True:
    lrange = float(input("Select the range of unmoving from 0 ~ 1.0 : "))
    if lrange >= 0 and lrange <= 1.0:
        break
    else:
        print("Error! Input the float 0 ~ 1.0 !")
        continue

# Set axis
while True: 
    axis = input("Select axis vertical to interface (x, y, z) : ")
    #axis = axis.split()
    if axis == "x" or axis == "X":
        axis_flag = 0
        break
    elif axis == "y" or axis == "Y":
        axis_flag = 1
        break
    elif axis == "z" or axis == "Z":
        axis_flag = 2
        break
    else:
        print("Error! Select from x, y and z !")
        continue


# read POSCAR
with open("POSCAR", "r") as f:
    lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].replace("\n", "")

flag_cartesian = True
if lines[7] == "Direct":
    flag_cartesian = False


# insert "Selectie dynamics"
lines.insert(7, "Selective dynamics")

# separate lattice & position
lattice_para = lines[0:9]
position = lines[9:]

for i in range(len(position)):
    position[i] = position[i].split()


# get z position and z length
lposition = []
for i in range(len(position)):
    lposition.append(float(position[i][axis_flag]))
lmax = max(lposition)
lmin = min(lposition)
length = lmax-lmin
print(length)
# set capability (TTT or FFF)
move = []
if flag_cartesian:
    for i in range(len(position)):
        if lposition[i] <= length*lrange+lmin:
            move.append("F F F") # don't move
        if lposition[i] > length*lrange+lmin:
            move.append("T T T") # can move

else:
     for i in range(len(position)):
        if lposition[i] <= lrange+lmin/length:
            move.append("F F F") # don't move
        if lposition[i] > lrange+lmin/length:
            move.append("T T T") # can move   

# write POSCAR_selected
with open("POSCAR_selected", "w") as f:
    for i in range(len(lattice_para)):
        f.write(lattice_para[i]+"\n")
    for i in range(len(position)):
        for j in range(3):
            f.write("{:.6f}\t".format(float(position[i][j])))
        f.write(move[i]+"\n")
