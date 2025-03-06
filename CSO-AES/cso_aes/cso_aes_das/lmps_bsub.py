import os
import sys


def lammps_bsub(dir, mlp_nums, mlp_encode_model,pwd):
    path = os.path.dirname(os.path.dirname(pwd))
    sys.path.append(os.path.join(os.path.dirname(__file__), os.path.join(path, 'init')))
    import bsub_script

    script_content_1 = f'''#!/bin/bash
{bsub_script.bsub_script_lmp_job_name}{dir}
'''
    script_content_1 = script_content_1 + bsub_script.bsub_script_lmp

    content_1 = '''pwd=$(pwd)
touch __start__
# Loop through directories in the current working directory
for dir in "$pwd"/*; do
    if [ -d "$dir" ]; then  # Check if it's a directory
        # Loop through subdirectories
        for sub_dir in "$dir"/*; do
            if [ -d "$sub_dir" ]; then  # Check if it's a directory
                echo "$sub_dir"
                cd "$sub_dir" 
                touch __start__
                # initial run 
    
                COMMAND="${COMMAND_0} -in lmp.in -var out_dump_file force.0.nc -var mlip_ini mlip_0.ini"
                $COMMAND > jobs_mlip_0.ini.out 2>&1

                # rerun
'''
    content_temp = ''

    if mlp_encode_model == False:
        for i in range(1, mlp_nums):
            content = f"""
                COMMAND="${{COMMAND_0}}  -in in_rerun.lmp -var out_dump_file force.{i}.nc  -var mlip_ini mlip_{i}.ini"
                $COMMAND > jobs_mlip_{i}.ini.out 2>&1
"""
            content_temp = content_temp + content
    content_temp = content_temp + """
                touch __ok__"""

    content_2 = """
            fi
        done
    fi
done
cd $pwd
touch __ok__
"""
    with open('bsub.lsf','w') as file:
        file.write(script_content_1+content_1+content_temp+content_2)


if __name__ =='__main__':
    #lammps_bsub('physics', 'CaTiO3', 'medium', '80', 3)
    lmp_exe = 'mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP /work/phy-huangj/app/interface-lammps-mlip-2/lmp_intel_cpu_intelmpi'
    lammps_bsub('qiming', 'CaTiO3', '38', '64', 3, lmp_exe,mlp_encode_model=False)
    #256G56c 56
