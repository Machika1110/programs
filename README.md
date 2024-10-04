# CALC_PROGRAM 説明
文責：内藤  
更新日：2024/07/30  
&emsp;CALC_PROGRAM内のフォルダの説明及びメモ。
## 0. アップロード方法
1. git pull origin main

2. git add .

3. git commit -m "コメント"

4. git push origin main

## 1. Structure
構造作成プログラム。出力はPOSCAR（VASP用）固定なので、**convert_file.py**で変換必要。

|プログラム名|結晶構造名|対応化合物|
---|---|---
SiO2_all.py|$\mathrm {SiO_2}$ の構造 | α-quartz, β-quartz, β-tridymite,  α-cristobalite, β-cristobalite, stishovite
bcc.py|体心立方構造|W, Cr, Fe
corundum.py|コランダム構造|$\mathrm {Al_2O_3}$
cubic.py|単体の立方晶系|Sb, Te
diamond.py|ダイヤモンド構造|$\mathrm {C, Si, Ge, Sn}$
diamond_mc.py|ダイヤモンド構造混晶|$\mathrm {SiC, SiGe, SiSn, GeSn}$
fcc.py|面心立方構造|Au, Ag, Cu, Fe, Al, Co
fluorite.py|蛍石型構造|$\mathrm {CaF_2, CeO_2}$
hcp.py|六方最密充填構造|Co, Ru
hexagonal.py|単体の六方晶系|$\mathrm {Co, Ru}$
perovskite.py|ペロブスカイト構造|$\mathrm {BaTiO_3, CaTiO_3}$
rocksalt.py|岩塩型構造|$\mathrm { Ge-Sb-Te, NaCl, MgO, CaO}$
sphalerite.py|閃亜鉛鉱型|$\mathrm {ZnS}$
rutile.py|ルチル型構造|$\mathrm {SiO_2, GeO_2}$
trigonal.py|単体の三方晶系|$\mathrm {Sb, Te}$
wurtzite.py|ウルツ型構造|$\mathrm {ZnO, ZnS, BeO, BN, GaN}$

### 使い方 (表の構造) 
1. 下の表を参考に、作成したい構造のプログラムを実行する。  

        $ python *****.py

2. 原子種を入力する。対応する数字を入力する。
3. サイズを入力する。例は以下の通り。

        Input number :  2 3 4

    x, y, z方向それぞれが2, 3, 4倍された超格子が作成できる。

4. POSCARファイルが作成される。
5. **convert_file.py**でファイル形式の変換を行う。

### 使い方（stack.py）
1. __同じファイル形式の__ くっつけたい構造を複数用意する。
2. 1つ目の構造を選択
   
        Please select 1st file :
3. 2つ目の構造を選択 
   
        Please select 2nd file : 
4. くっつける方向をx, y, z から指定。**小文字で入力**
     
        Please input direction to stack : 
5. 2つの構造間の距離を決定。単位はオングストローム。
   
        Please select distance among two structure : 
6. くっつけた後のファイル名を入力。
   
        Please input the name of the converted file:

### cycro-SiGe
Si：Ge＝1：1の混晶のプログラム。SiGe混晶の構造を作成後、DFSアルゴリズムに基づき、SiおよびGeの**単体の六員環構造**を含んでいるかどうかも出力。

## 2. Tool_Allegro
### calc_temp
&emsp;AllegroでLAMMPS計算を行う際のテンプレート

### comb_train.py
&emsp;訓練データの拡張.xyzファイル（vasprun.xmlをxml_to_extxyz.pyで変換）をまとめて一つのファイルへ格納する。

## 3. Tool_HDNNP
### comb_inputdata.py
&emsp;HDNNP用の訓練データinput.data形式のファイルを結合する。

## 4. Tool_LAMMPS
### trace_lmp.py
log.lammpsで出力された計算過程をグラフ化。  
表はExcel、グラフは.jpg形式でフォルダにまとめて出力。

### ratio-gaussian.py
calc.finalなどのlammpsのdumpファイルから、ｚ方向の各原子の濃度分布を縦軸、ｚ座標を横軸に取ったグラフ（.png）と表（.xlsx）を出力する。

## 5. Tool_VASP
### calc_mlff_***.py
ctifor : 閾値を継続  
 shot  : 原子を追加し、発射する  
 無印  ：閾値を毎回初期値に

### calc_scale.py
プログラム内の値に応じ、原子の構造を等方的にスケーリングして計算する。  
EV_graph.pyと使うことで、EV図を作成できる。  

### EV_graph.py
計算結果の体積とエネルギーを取得し、EV図を作成できる。  
calc_scale.pyと組み合わせでEV図を作成できる。

### kpoints.py
POSCARを読み、適切なKPOINTSを生成。格子が最小単位でない場合は警告文が出る。

### plot_band.py
vasprun.xmlを読み、電子のバンド図を作成。なお、CHGCARを読んだ計算（ICHARGE=-11）にのみ限定

### PVTE_graph.py
VASPの計算結果から、圧力、体積、温度、エネルギーのグラフを作成。

### SD_single.py / SD_double.py
POSCAR内の原子を一部固定する。singleはz方向の小さい方、doubleはz方向の上下を固定する。

### xml_to_extxyz.py
第一原理MDのvasprun.xmlから、Allegro・NequIP用の訓練データファイルを作成。freqで取り出す間隔を決められる。出力ファイル名の数字は訓練データの数。

### xml_to_inputdata.py
第一原理MDのvasprun.xmlから、HDNNP用の訓練データファイルを作成。freqで取り出す間隔を決められる。出力ファイル名の数字は訓練データの数。

## 6. other
### convert-data.py  
n2p2用の訓練データファイルinput.dataとAllegro・NequIP用訓練データファイル拡張xyzファイルを変換するプログラム。

### extxyz-to-input.py  
前述の通り

### input-to-extxyz.py  
前述の通り