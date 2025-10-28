from ase.io import read, write
import os
import shutil
import yaml
import sys
import glob

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

def read_poscar(poscar_file="POSCAR"):
    """读取POSCAR文件，返回元素列表"""
    try:
        with open(poscar_file, 'r') as f:
            lines = f.readlines()

        # 找到元素行（通常是第6行）
        elements_line = None
        for i in range(3, min(10, len(lines))):
            line = lines[i].strip()
            # 简单的元素检测：包含大写字母且不全是数字
            if line and any(c.isupper() for c in line) and not all(c.isdigit() or c.isspace() for c in line):
                elements_line = line
                break

        if not elements_line:
            raise ValueError("Can't find the element row in the POSCAR file!")

        elements = elements_line.split()
        #print(f"在POSCAR中找到元素: {elements}")
        return elements

    except Exception as e:
        raise ValueError(f"There is an error reading POSCAR: {e}!")



def find_potcar_file(element, potcar_dir):
    """使用模糊搜索查找元素的POTCAR文件"""
    # 可能的搜索模式
    search_patterns = [
        f"{element}*/POTCAR",  # 元素文件夹下的POTCAR文件
        f"POTCAR_{element}",  # 直接的POTCAR_元素文件
        f"*/POTCAR_{element}",  # 任何子目录下的POTCAR_元素文件
        f"{element}*",  # 以元素名开头的文件（备用）
        f"*{element}*",  # 包含元素名的文件（备用）
    ]

    for pattern in search_patterns:
        search_path = os.path.join(potcar_dir, pattern)
        matches = glob.glob(search_path)

        # 过滤出文件（排除目录）
        file_matches = [m for m in matches if os.path.isfile(m)]

        if file_matches:
            # 优先选择更精确的匹配
            for match in file_matches:
                filename = os.path.basename(match)
                parent_dir = os.path.basename(os.path.dirname(match))

                # 最优先：元素文件夹下的POTCAR文件
                if filename == "POTCAR" and parent_dir.startswith(element):
                    return match
                # 其次：直接的POTCAR_元素文件
                elif filename == f"POTCAR_{element}":
                    return match

            # 返回第一个匹配的文件
            return file_matches[0]

    return None


def create_potcar(poscar_file, potcar_dir, output_file="POTCAR"):
    """合成POTCAR文件"""

    if not os.path.exists(potcar_dir):
        raise ValueError(f"error! pseudo potential ('{potcar_dir}') isn't exist !")

    elements = read_poscar(poscar_file)
    total_size = 0
    found_files = []

    with open(output_file, 'wb') as outfile:
        for element in elements:
            # 使用模糊搜索查找元素的POTCAR文件
            potcar_file = find_potcar_file(element, potcar_dir)

            if potcar_file and os.path.isfile(potcar_file):
                filename = os.path.basename(potcar_file)
                #print(f"添加 {element} 的POTCAR: {filename}")
                try:
                    with open(potcar_file, 'rb') as infile:
                        content = infile.read()
                        outfile.write(content)
                        total_size += len(content)
                    found_files.append((element, filename))
                except Exception as e:
                    print(f"There is an error reading {element} {potcar_file}.  {e}!")
            else:
                raise ValueError(f"error!: Can't find {element} POTCAR file in {potcar_dir}")

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
        r2scan_path = os.path.join(path, 'init', 'cp2k_r2scan.inp')
        if os.path.exists(r2scan_path):
            shutil.copy(r2scan_path, destination_folder)
        chek_file_exist(INPUT_path, scf_cal_engine)
        shutil.copy(INPUT_path, destination_folder)
    if scf_cal_engine == 'vasp':
        INPUT_path = os.path.join(path, 'init', 'INCAR')
        KPOINTS_path = os.path.join(path, 'init', 'KPOINTS')
        if os.path.exists(KPOINTS_path):
            shutil.copy(KPOINTS_path, destination_folder)
        chek_file_exist(INPUT_path, scf_cal_engine)
        shutil.copy(INPUT_path, destination_folder)




