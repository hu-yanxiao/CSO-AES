from ase.io import read, write
import os
import shutil
import yaml

def find_vasp_files_in_current_directory():
    current_directory = os.getcwd()  # 获取当前目录
    vasp_files = [file for file in os.listdir(current_directory) if file.endswith('.vasp')]
    return vasp_files
def get_upper_directory(level, current_path):
    for _ in range(level):
        current_path = os.path.dirname(os.path.abspath(current_path))
    return current_path

def poscar2STRU(dir):
    if not os.path.exists('POSCAR'):
        poscar = find_vasp_files_in_current_directory()[0]
        shutil.copy(poscar,'POSCAR')
    path = get_upper_directory(5, dir)
    pp_orb_dic = os.path.join(path, 'init', 'pp_orb_dic.yaml')


    with open(pp_orb_dic, 'rb') as file:
        yaml_data = yaml.safe_load(file)
        pp_dic = yaml_data['pp_dic']
        orb_dic = yaml_data['orb_dic']

    pwd = os.getcwd()
    cs_vasp = os.path.join(pwd, "POSCAR")
    #cs_vasp = os.path.join(pwd, "POSCAR")
    cs_atoms = read(cs_vasp, format="vasp")
    cs_stru = os.path.join(pwd, "STRU")

    atom_symbols = list(set(cs_atoms.get_chemical_symbols()))

    # # 需要设置原子磁矩时，可添加代码：
    # cs_atoms.set_initial_magnetic_moments([1.0,1.0,1.0,2.0])  # 后面列表里列出每个原子的磁矩；

    # # 如果需要设置每个原子每个方向的磁矩，则需要设置为二维数组：
    # cs_atoms.set_initial_magnetic_moments([[1.0,1.0,1.0],
    #                                        [1.0,1.0,2.0],
    #                                        [1.0,1.0,3.0],
    #                                        [1.0,1.0,4.0]])
    pp = {}
    basis = {}
    for atom in atom_symbols:
        pp[atom] = pp_dic[atom]
        basis[atom] = orb_dic[atom]

    write(cs_stru, cs_atoms, format="abacus", pp=pp, basis=basis)

# def modify_KPT(input, out_name):
#     with open(input, 'r') as input_file, open(out_name, 'w') as output_file:
#         lines = input_file.readlines()
#         lines[0] = 'K_POINTS\n'
#         lines[3] = lines[3].strip() + '   '+ lines[4]
#         for line in lines[:4]:
#             output_file.write(line)
#
# os.system('(echo 102; echo 2; echo 0.03)|vaspkit')
# modify_KPT('KPOINTS',"KPT")
# os.system('rm INCAR POTCAR KPOINTS ase_sort.dat')

def chek_file_exist(file_path,scf_cal_engine):
    if not os.path.exists(file_path):
        raise ValueError(f'{scf_cal_engine}:{file_path} is not exist!')

def INPUT(dir,destination_folder,scf_cal_engine):
    path = get_upper_directory(5, dir)
    if scf_cal_engine == 'abacus':
        INPUT_path = os.path.join(path,'init','INPUT')
        chek_file_exist(INPUT_path, scf_cal_engine)
        shutil.copy(INPUT_path,destination_folder)
    if scf_cal_engine == 'cp2k':
        INPUT_path = os.path.join(path, 'init', 'cp2k.inp')
        chek_file_exist(INPUT_path, scf_cal_engine)
        shutil.copy(INPUT_path, destination_folder)
    if scf_cal_engine == 'vasp':
        INPUT_path = os.path.join(path, 'init', 'INCAR')
        chek_file_exist(INPUT_path, scf_cal_engine)
        shutil.copy(INPUT_path, destination_folder)




