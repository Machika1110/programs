import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.core.periodic_table import Element
import shutil
import math

# standard of big error of energy / atom [eV]
standard_error = 0.5

# テストデータの割合
ratio_test = 0.1

# 現在のディレクトリを取得
current_directory = os.getcwd()


# xyzファイル名
input_xyz = 'train-seq.xyz'

# LAMMPSの実行ファイル
lammps_executable = 'lmp_7net'

# 出力用ディレクトリ名
output_dir = './structures'

# 出力ディレクトリの作成
os.makedirs(output_dir, exist_ok=True)

# 格子定数の変換
def convert_basis(basis):
    basis = np.array([[float(basis[0][0]), float(basis[0][1]), float(basis[0][2])],
                      [float(basis[1][0]), float(basis[1][1]), float(basis[1][2])],
                      [float(basis[2][0]), float(basis[2][1]), float(basis[2][2])]])

    num_basis, dim = basis.shape

    lengths = np.linalg.norm(basis, axis=1)

    unit_basis = basis / lengths[:, np.newaxis]

    angles = np.zeros((num_basis, num_basis))
    for i in range(num_basis):
        for j in range(i + 1, num_basis):
            angle = np.arccos(np.dot(unit_basis[i], unit_basis[j]))
            angles[i][j] = angle
            angles[j][i] = angle

    A, B, C = lengths[0], lengths[1], lengths[2]
    alpha ,beta, gamma = angles[1][2], angles[0][2], angles[0][1]

    ax = A
    bx = B * math.cos( gamma )
    by = B * math.sin( gamma )
    cx = C * math.cos( beta )
    cy = ( B * C * math.cos( alpha ) - bx * cx ) / by
    cz = np.sqrt( C ** 2 - cx ** 2 - cy ** 2 )
    converted_basis =[[ax, bx, cx],
                      [0, by, cy],
                      [0, 0, cz]]

    return converted_basis

# 座標変換
def convert_coor(lattice, coordinate, changed_lattice):
    coor = [] 
    changed_basis = np.array([[changed_lattice[0][0], changed_lattice[0][1], changed_lattice[0][2]], 
                              [changed_lattice[1][0], changed_lattice[1][1], changed_lattice[1][2]], 
                              [changed_lattice[2][0], changed_lattice[2][1], changed_lattice[2][2]]])
    A = np.array([float(lattice[0][0]), float(lattice[0][1]), float(lattice[0][2])])
    B = np.array([float(lattice[1][0]), float(lattice[1][1]), float(lattice[1][2])])
    C = np.array([float(lattice[2][0]), float(lattice[2][1]), float(lattice[2][2])])
    
    Volume = abs(np.dot(np.cross(A, B), C))
    
    tmp_array = np.array([[np.cross(B, C)[0], np.cross(B, C)[1], np.cross(B, C)[2]], 
                          [np.cross(C, A)[0], np.cross(C, A)[1], np.cross(C, A)[2]],
                          [np.cross(A, B)[0] , np.cross(A, B)[1], np.cross(A, B) [2]]])
    tmp_array = tmp_array / Volume

    for i in range(len(coordinate)):
        coor_before = np.array([[float(coordinate[i][0])], [float(coordinate[i][1])], [float(coordinate[i][2])]])
        coor_after = np.dot(np.dot(changed_basis, tmp_array), coor_before)
        tmp_coor = [coor_after[0][0], coor_after[1][0], coor_after[2][0]]
        coor.append(tmp_coor)
    
    return coor

# XYZファイルの読み込み
def read_xyz(file_path):

    # ファイルを読み込みdataに格納
    with open(file_path, 'r') as file:
        data = file.readlines()
    
    # 構造のデータを格納するための配列
    structures = []
    index = 0

    while index < len(data):
        # 原子数をnum_atomsに格納
        num_atoms = int(data[index].strip())

        # 格子やエネルギーが書かれている行をcommentに格納
        comment = data[index + 1].strip()

        # 原子座標と原子に働く力をatomsに格納
        atoms = data[index + 2:index + 2 + num_atoms]

        # 格子の情報をcommentから取り出す
        lattice = comment.split()[:9]

        # 各情報をstructuresに格納
        structures.append((num_atoms, comment, atoms, lattice))

        # 次の構造を読み込むため、indexの値を変更
        index += 2 + num_atoms
    return structures

# LAMMPSの構造ファイル作成
def write_lammps_structure(file_path, structure):
    num_atoms, comment, atoms, lattice = structure

    # 構造ファイルを作成
    with open(file_path, 'w') as file:
        file.write("LAMMPS data file\n\n")
        file.write(f"{num_atoms} atoms\n")

        # 原子種を読み込みatom_typesに格納
        atom_types = []
        for i, atom in enumerate(atoms):
            parts = atom.split()
            atom_types.append(parts[0])
        
        # 原子種を取得
        atom_type = sorted(set(atom_types), key=atom_types.index)
        file.write(str(len(atom_type)) + " atom types\n")
        
        # 格子ベクトルの大きさを書き込む
        lattice[0] = lattice[0].strip("Lattice=\"")
        lattice[-1] = lattice[-1].strip("\"")

        # 格子ベクトルを二次元配列に格納
        new_lattice = [[lattice[0], lattice[1], lattice[2]], [lattice[3], lattice[4], lattice[5]], [lattice[6], lattice[7], lattice[8]]]

        # 対角成分が0出ない場合flagを立てる
        for i in range(len(new_lattice)):
            for j in range(len(new_lattice[i])):
                if i != j and float(new_lattice[i][j]) != 0.0:
                    flag = 1
                    break

        # 座標を二次元配列に格納
        coordinate = []
        for i in range(len(atoms)):
            tmp = atoms[i].split()
            coordinate.append([tmp[1], tmp[2], tmp[3]])

        # 格子の形によって書式を変える
        if flag == 0:
            file.write("0.0 " + lattice[0] + " xlo xhi\n")
            file.write("0.0 "+ lattice[4] + " ylo yhi\n")
            file.write("0.0 " + lattice[-1] + " zlo zhi\n\n")
            converted_coor = coordinate
        elif flag == 1:
            converted_basis = convert_basis(new_lattice)
            converted_coor = convert_coor(new_lattice, coordinate, converted_basis)
            xhi = converted_basis[0][0]
            yhi = converted_basis[1][1]
            zhi = converted_basis[2][2]
            xy = converted_basis[0][1]
            xz = converted_basis[0][2]
            yz = converted_basis[1][2]
            file.write(" 0.00 {} xlo xhi\n".format(xhi))
            file.write(" 0.00 {} ylo yhi\n".format(yhi))
            file.write(" 0.00 {} zlo zhi\n".format(zhi))
            file.write(" {} {} {} xy xz yz\n\n".format(xy ,xz ,yz))

        file.write("Masses\n\n")

        # 対応する原子種の質量を取得
        for i, j in enumerate(atom_type):
            atom_data = Element(j)
            atom_mass = str(atom_data.atomic_mass)
            file.write(str(i+1) + " " + str(atom_mass.split()[0]) + "\n")
        file.write("\n")

        # 原子の座標を書き込む
        file.write("Atoms\n\n")
        for i in range(len(converted_coor)):
            # 通し番号を書き込む
            file.write(f"{i + 1} ")

            #各原子固有の番号を取得
            for p, q in enumerate(atom_type):
                if atom_types[i] == q: 
                    tmp = p
                    break

            file.write(f"{tmp+1} {converted_coor[i][0]} {converted_coor[i][1]} {converted_coor[i][2]}\n")
    
    # in.reg-7netを編集
    # in.reg-7net_newを作成
    shutil.copy("in.reg-7net", "in.reg-7net-new")

    # in.reg-7net_newの中身を編集
    # in.reg-7net_newを読み込む
    with open("in.reg-7net-new", "r") as f:
        lines = f.read()

    # 構造ファイル名を読み込ませるためにread_dataを編集
    lines = lines.replace("read_data    tmp.lmp", "read_data    " + file_path)
    tmp = ""

    # ex) atom_type = [Ni, Ge)のとき、tmp = "Ge Ni"
    for i in atom_type:
        tmp += " " + i
        
    print(tmp)

    # 計算対象の系に含まれる元素に合わせて、pair_coeffを編集
    lines = lines.replace("pair_coeff", "pair_coeff   * * deployed_serial.pt {}".format(tmp))

    # in.reg-7net_newへ書き込む
    with open("in.reg-7net-new", "w") as f:
        f.write(lines)

# 1原子あたりのエネルギー計算
def energy_per_atom(energy, num_atoms):
    return energy / num_atoms

# LAMMPSの計算を実行
def run_lammps(log_file):
    lammps_input_template = "in.reg-7net-new"
    subprocess.run([lammps_executable, '-in', lammps_input_template, '-log', log_file])

### メイン処理 ###
structures = read_xyz(input_xyz)
dft_energies = []
lammps_energies = []

for i, structure in enumerate(structures):
    num_atoms, comment, atoms, lattices = structure
    
    # DFTエネルギーを抽出（stress以降を無視）
    energy_part = comment.split('stress=')[0].strip()

    # energy=●●● の形式のため、energy=を削除
    dft_energy = float(energy_part.split()[10].strip("energy="))

    # 1原子当たりのエネルギーを取得しdft_energiesに格納
    dft_energy_per_atom = energy_per_atom(dft_energy, num_atoms)
    dft_energies.append(dft_energy_per_atom)
    
    # 作成する構造ファイルの名前を決定
    structure_file = os.path.join(output_dir, f'structure_{i}.data')

    # 構造ファイルの作成
    write_lammps_structure(structure_file, structure)
    
    # ログファイルの名前を決定
    log_file = os.path.join(output_dir, f'log_{i}.lammps')

    # Lammpsの実行
    run_lammps(log_file)
    
    # logファイルからポテンシャルエネルギーを取得
    with open(log_file, 'r') as log_file:
        lines = log_file.readlines()
        flag = 0
        for j in range(len(lines)):
            # Per MPI rank memory　と書かれた行の二行後にポテンシャルエネルギーが記載される
            if "Per MPI rank memory" in lines[j]:
                e_pot_index = lines[j+1].split().index("PotEng")
                lammps_energy = float(lines[j+2].split()[e_pot_index])
                lammps_energy_per_atom = energy_per_atom(lammps_energy, num_atoms)
                lammps_energies.append(lammps_energy_per_atom)     
                flag = 1
                break
        
        # 計算が動かない等不具合が起きた時用
        if flag == 0:
            lammps_energies.append(0.0)




# 回帰グラフの作成（訓練データ）
plt.figure(figsize=(10, 5))
plt.scatter(dft_energies, lammps_energies, label='Data Points')
plt.plot(dft_energies, dft_energies, color='red', linestyle='--', label='Ideal Fit')
plt.xlabel('DFT Energy per Atom (eV)')
plt.ylabel('LAMMPS Energy per Atom (eV)')
plt.title('Regression Plot of DFT vs. LAMMPS Energy of train data')
plt.legend()
plt.grid(True)
plt.show()
plt.savefig("reg.png")



# conserve to .txt file
with open("regression.txt", "w") as f:
    f.write("energy per atom [eV] of train data\n")
    f.write("DFT energy  LAMMPS energy\n")
    for i in range(len(dft_energies)):
        f.write("{} {}\n".format(dft_energies[i], lammps_energies[i]))

# conserve to .txt file of datum which show big error based on "standard_error"
with open("dangerous_data.txt", "w") as f:
    f.write("number of conf.   abs error\n")
    for i in range(len(dft_energies)):
        if abs(float(dft_energies[i])-float(lammps_energies[i])) > standard_error:
            f.write("{}   {}\n".format(i+1, abs(float(dft_energies[i])-float(lammps_energies[i]))))
