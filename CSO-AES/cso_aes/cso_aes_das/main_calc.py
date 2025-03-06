from ase.io import iread, write
import os
import shutil
from .gen_calc_file import poscar2STRU,INPUT
from tqdm import trange
import yaml
import sys

def check_and_modify_calc_dir(pwd,xyz_list,calc_dir_num):
    parameter_yaml = os.path.join(pwd,'parameter.yaml')
    end_yaml = os.path.join(pwd,'end.yaml')
    xyz_num = len(xyz_list)

    if xyz_num < calc_dir_num:
        new_calc_dir_num = xyz_num
        with open(parameter_yaml, 'r') as file:
            data_1 = yaml.safe_load(file)
        dft = data_1['dft']
        dft['calc_dir_num'] = new_calc_dir_num
        with open(parameter_yaml, 'w') as file:
            yaml.safe_dump(data_1, file, default_flow_style=False)

        with open(end_yaml, 'r') as file:
            data_2 = yaml.safe_load(file)
        dft = data_2['dft']
        dft['calc_dir_num'] = new_calc_dir_num
        with open(end_yaml, 'w') as file:
            yaml.safe_dump(data_2, file, default_flow_style=False)
        return new_calc_dir_num
    else:
        return calc_dir_num

def mkdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

####数据集划分#######
def mkdir_Placefile(dir,num_posacr,calc_dir_num,xyz,scf_cal_engine):
    quotient = num_posacr // calc_dir_num
    remainder = num_posacr % calc_dir_num
    temp = [quotient] * calc_dir_num
    divide = [a + 1 if index < remainder else a for index, a in enumerate(temp)]
    range_list = [[sum(divide[:a + 1]) - divide[a] + 1, sum(divide[:a + 1])] for a in range(len(divide))]
    for i in trange(1,calc_dir_num+1):
        sub_dir = os.path.join(dir, 'dir_' + str(i))
        mkdir(sub_dir)
        index_poscar = [a for a in range(range_list[i - 1][0], range_list[i - 1][1] + 1)]

        for ii in index_poscar:
            sub_sub_dir = os.path.join(sub_dir, str(ii))
            mkdir(sub_sub_dir)
            if scf_cal_engine == 'abacus':
                write(os.path.join(sub_sub_dir,'POSCAR'), xyz[ii-1], format='vasp')
                os.chdir(sub_sub_dir)
                poscar2STRU(dir)
                INPUT(dir,sub_sub_dir,scf_cal_engine)
                os.remove('ase_sort.dat')
            elif scf_cal_engine == 'cp2k':
                atoms = xyz[ii - 1]
                write(os.path.join(sub_sub_dir, 'cp2k.xyz'), atoms, format='extxyz')
                llist = atoms.cell.cellpar()
                os.chdir(sub_sub_dir)
                INPUT(dir, sub_sub_dir, scf_cal_engine)
                cell = f'     ABC {llist[0]}  {llist[1]}  {llist[2]} \n'
                angles = f'     ALPHA_BETA_GAMMA {llist[3]}  {llist[4]}  {llist[5]}\n'
                cp2k_inp = 'cp2k.inp'
                with open(cp2k_inp, 'r') as f:
                    lines = f.readlines()
                    for index, line in enumerate(lines):
                        tt = index
                        if '&CELL' in line:
                            break

                    lines[tt + 1] = cell
                    lines[tt + 2] = angles
                with open(cp2k_inp, 'w') as file:
                    file.writelines(lines)
            elif scf_cal_engine == 'vasp':
                b = xyz[ii-1]
                sorted_indices = sorted(range(len(b)), key=lambda i: b[i].symbol)
                sorted_atoms = b.__class__()
                sorted_atoms.cell = b.cell
                sorted_atoms.pbc = b.pbc
                for i in sorted_indices:
                    sorted_atoms.append(b[i])
                write(os.path.join(sub_sub_dir,'POSCAR'), sorted_atoms, format='vasp')
                os.chdir(sub_sub_dir)
                INPUT(dir,sub_sub_dir,scf_cal_engine)
                os.system(f'(echo 103)|vaspkit')
            else:
                raise ValueError(f'{scf_cal_engine} is not exist!')

def main_calc(atom_list,calc_dir_num,path_main):
    #file = [f for f in os.listdir() if f.endswith('.xyz')][0]
    # calc_dir_num=13 #计算文件数目
    #num = 1  #1:scf 2:cell-relax
    #server = 2 #1:qiming 2:wulixi

    path = os.path.dirname(os.path.dirname(path_main))
    jobname = os.path.basename(path)
    sys.path.append(os.path.join(os.path.dirname(__file__), os.path.join(path, 'init')))
    import bsub_script

    pwd = os.getcwd()
    dir = os.path.join(pwd,'filter')
    mkdir(dir)

    num_posacr = len(atom_list)

    mkdir_Placefile(dir,num_posacr,calc_dir_num,atom_list,bsub_script.scf_cal_engine)

    os.chdir(pwd)
    jobs_script = os.path.join(pwd,'jobs_script')
    mkdir(jobs_script)

    quotient = num_posacr // calc_dir_num
    remainder = num_posacr % calc_dir_num
    temp = [quotient] * calc_dir_num
    divide = [a + 1 if index < remainder else a for index, a in enumerate(temp)]
    range_list = [[sum(divide[:a + 1]) - divide[a] + 1, sum(divide[:a + 1])] for a in range(len(divide))]

    ####创建任务提交脚本#####
    for i in range(1,calc_dir_num+1):
        range1,range2 = range_list[i - 1][0], range_list[i - 1][1]
        #print(range1,range2)
        content_1 = f'''#!/bin/bash
{bsub_script.bsub_script_scf_job_name} {jobname}_{str(i)} 
'''
        content_1 = content_1 + bsub_script.bsub_script_scf

        content_2 = f'''dir_1="filter"
dir_2="{'dir_' + str(i)}"
cd ../$dir_1/$dir_2

path=$(pwd)

for item in {{{range1}..{range2}}}; do
    cd $path/$item
    start_time=$(date +%s.%N)
    touch __start__
    $COMMAND_std > logout 2>&1
    touch __ok__
    
    end_time=$(date +%s.%N)
    runtime=$(echo "$end_time - $start_time" | bc)
    cd $path
    cd ..
    echo "job_{str(i)}_$item total_runtime:$runtime s" >> time.txt
done'''
        with open(os.path.join(jobs_script,f'bsub_{i}.lsf'), 'w') as file:
            file.write(content_1+content_2)
        os.chdir(jobs_script)
        #os.system(f'dos2unix bsub_{i}.lsf')

    os.chdir(pwd)

    start_calc = f'''import os
start = 1
end = {calc_dir_num}

pwd = os.getcwd()
jobs_script = os.path.join(pwd,'jobs_script')
os.chdir(jobs_script)
for i in range(start,end+1):
    #print(f'python bsub_{{i}}.lsf')
    if not os.path.exists(f'bsub_{{i}}.lsf'):
        print('errors')
    os.system(f'{bsub_script.start_calc_command} bsub_{{i}}.lsf')'''

    with open('start_calc.py', 'w') as f:
        f.write(start_calc)
    return num_posacr