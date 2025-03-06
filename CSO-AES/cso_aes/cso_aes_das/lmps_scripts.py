from ase.io import read, write
import os
import shutil
from ase.data import atomic_numbers
from ase.data import atomic_masses
from ase.data import chemical_symbols
import sys
import random
import importlib.util
import sys

def main_pos2lmp(ele_, size, ele_model):
    if os.path.exists('POSCAR'):
        POSCAR = 'POSCAR'
    else:
        POSCAR = [f for f in os.listdir() if f.endswith('.vasp')][0]
    a = read(POSCAR)
    a = a.repeat(size)

    if ele_model == 1:
        ele = [atomic_numbers[i] for i in ele_]
        ele = sorted(ele)
        specorder = [chemical_symbols[i] for i in ele]

    elif ele_model == 2:
        specorder = ele_
    else:
        raise ValueError('ele_model must be 1 or 2!')

    write('data.in', a, format='lammps-data', masses=True, specorder=specorder, force_skew=True)

def lammps_scripts(ensemble,temp,mlp_nums,mlp_encode_model,pwd):
    path = os.path.dirname(os.path.dirname(pwd))
    sys.path.append(os.path.join(os.path.dirname(__file__), os.path.join(path, 'init')))
    import lmp_in
    random_number = random.randint(1000, 1000000)
    for i in range(0,mlp_nums):
        content = f'''mtp-filename current_{i}.mtp
select FALSE'''

        with open(f'mlip_{i}.ini','w') as file:
            file.write(content)

    pre_data_in = f"""variable T equal {temp}
    variable random equal {random_number}
    """

    if ensemble=='npt':
        lmp_in__ = pre_data_in + lmp_in.npt
    elif ensemble=='nvt':
        lmp_in__ = pre_data_in + lmp_in.nvt
    else:
        fix = f'{ensemble}: error'
        raise ValueError(fix)

    return_lmp = f"""variable input_dump_file string "force.0.dump"
# variable out_dump_file string 
# variable mlip_ini string 

units metal
boundary p p p 
atom_style atomic

box tilt large
read_data data.in

pair_style mlip ${{mlip_ini}}
pair_coeff * *

neighbor 2.0 bin 
neigh_modify delay 10 check yes 

reset_timestep 0

compute pe all pe/atom
dump 1 all custom 1 ${{out_dump_file}} id type x y z fx fy fz c_pe
dump_modify 1 sort id

rerun ${{input_dump_file}} dump x y z"""

    with open('lmp.in', 'w') as file:
        file.write(lmp_in__)
    if mlp_encode_model == False:
        with open('in_rerun.lmp', 'w') as file:
            file.write(return_lmp)




