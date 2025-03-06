import os.path
from ase.io import iread,write
from .work_dir import bsub_dir
import glob
import pickle

def mkdir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

'''生成scf相关目录，记录为收敛的结构(stru.pkl)，返回需要scf的结构'''
def scf_lammps_data(delete_filter_select,pwd):
    total_atom_list = []
    scf_lammps_data_path = os.path.join(pwd,'scf_lammps_data')
    scf = os.path.join(scf_lammps_data_path,'scf')
    #total_sample_filter_path = os.path.join(scf,'total_sample_filter.xyz')

    mkdir(scf_lammps_data_path)
    mkdir(scf)
    stru_list = []
    stru_pkl = os.path.join(pwd,'stru.pkl')
    for dir in bsub_dir(pwd):
        os.chdir(dir)
        delete = glob.glob('*filter*')
        sample = [i for i in delete if 'sample' in i]
        delete = [i for i in delete if 'remain' in i]
        if len(sample) != 0:
            stru_list.append(os.path.basename(dir)+'.vasp')
            total_atom_list = total_atom_list + list(iread(sample[0]))
        if delete_filter_select == 'yes':
            for i in delete:
                if os.path.exists(i):
                    os.remove(i)
    with open(stru_pkl, 'wb') as f:
        pickle.dump(stru_list, f)
    #write(total_sample_filter_path,total_atom_list,format='extxyz')
    return total_atom_list


