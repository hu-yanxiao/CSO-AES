'''
Ref
taiyi: mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP /work/phy-huangj/app/interface-lammps-mlip-2/lmp_intel_cpu_intelmpi
physics: mpirun /share/home/wuyb/appbins/interface-lammps-mlip-2/lmp_mpi
qiming: mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP /work/phy-huangj/app/interface-lammps-mlip-2/lmp_intel_cpu_intelmpi


#SBATCH --job-name lmp_{dir}
#SBATCH --partition {queue}
##SBATCH --nodelist c0[01-40]
#SBATCH --ntasks {cores}  #number of core
#SBATCH --qos=840cpu

ulimit -s unlimited
ulimit -l unlimited
COMMAND_0="{lmp_exe}
'''

'''
str = f"""#!/bin/bash
#SBATCH --job-name lmp_{current_num}
#SBATCH --partition {train_queue}
##SBATCH --nodelist c0[01-40]
#SBATCH --ntasks {train_cores}  #number of core
#SBATCH --qos=840cpu
module load mkl/2022.1.0 mpi/2021.6.0 
ulimit -s unlimited
ulimit -l unlimited
'''



train_sus_queue = 33
train_sus_cores = 40
train_sus_ptile = 40

'''train_sus_exe'''
sus2_mlp_exe = '/work/phy-huangj/app/SUS2-MLIP/bin/mlp-sus2'
# original_COMMAND = f"mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP {sus2_mlp_exe} train hyx.mtp ../train.cfg  --stdd-weight=0.0   --std-weight=0.0000   --max-iter=2000  --curr-pot-name="
# subsequent_COMMAND = f"mpiexec.hydra -machinefile $LSB_DJOB_HOSTFILE -np $NP {sus2_mlp_exe} train hyx.mtp ../train.cfg   --init-params=same --do-samp=false  --stdd-weight=0.0   --std-weight=0.0000   --max-iter=1000  --curr-pot-name="
original_COMMAND = f"mpirun {sus2_mlp_exe} train hyx.mtp ../train.cfg  --stdd-weight=0.0   --std-weight=0.0000   --max-iter=2000  --curr-pot-name="
subsequent_COMMAND = f"mpirun {sus2_mlp_exe} train hyx.mtp ../train.cfg   --init-params=same --do-samp=false  --stdd-weight=0.0   --std-weight=0.0000   --max-iter=1000  --curr-pot-name="


lmp_queue = 33
lmp_cores = 40
lmp_ptile = 40
lmp_exe = 'mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP /work/phy-huangj/app/interface-lammps-mlip-2/lmp_intel_cpu_intelmpi'

scf_queue = 33
scf_cores = 40
scf_ptile = 40
scf_cal_engine = 'abacus'       # abacus and cp2k

##BSUB -J lmp_
##BSUB -J lmp_
##BSUB -J

##SBATCH --job-name lmp_
##SBATCH --job-name

bsub_script_train_sus_job_name = f'''
#BSUB -J lmp_'''

bsub_script_lmp_job_name = f'''
#BSUB -J lmp_'''

bsub_script_scf_job_name = f'''
#BSUB -J'''

bsub_script_train_sus = f'''
#BSUB -q {train_sus_queue}
#BSUB -n {train_sus_cores}
#BSUB -R "span[ptile={train_sus_ptile}]"
#BSUB -e %J.err
#BSUB -o %J.out
hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
module load mpi/2021.6.0 compiler/2022.1.0 mkl/2022.2.0  
'''

bsub_script_lmp = f'''
#BSUB -q {lmp_queue}
#BSUB -n {lmp_cores}
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile={lmp_ptile}]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
COMMAND_0="{lmp_exe}"
'''


if scf_cal_engine == 'abacus':
    bsub_script_scf = f'''
#BSUB -q {scf_queue}
#BSUB -n {scf_cores}
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile={scf_ptile}]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
module load compiler/2022.1.0 mpi/2021.6.0 mkl/2022.2.0
export PATH=/work/phy-huangj/apps/il/abacus/3.6.5/bin:$PATH
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP abacus"
'''

#     bsub_script_scf = f'''
# #SBATCH --partition {scf_queue}
# ##SBATCH --nodelist c0[01-40]
# #SBATCH --ntasks {scf_cores}  #number of core
# #SBATCH --qos=840cpu
# #####
# module load mkl/2022.1.0 mpi/2021.6.0 compiler/2022.1.0
# #export LD_LIBRARY_PATH=/share/home/wuyb/opts/apps/hdf5-1.10.9/lib:$LD_LIBRARY_PATH
# conda activate abacus_env
# ulimit -s unlimited
# ulimit -l unlimited
# export PATH=/share/home/wuyb/appbins/abacus-3.6.5/bin:$PATH
# COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP abacus"
# '''


elif scf_cal_engine == 'cp2k':
    bsub_script_scf = f'''
#BSUB -q {scf_queue}
#BSUB -n {scf_cores}
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile={scf_ptile}]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
module purge
module load cp2k/2024.1_oneapi-e5
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP cp2k.popt -i cp2k.inp"
'''

###vasp 生成POTCAR文件采用的是vaspkit软件
elif scf_cal_engine == 'vasp':
    bsub_script_scf = f'''
#BSUB -q {scf_queue}
#BSUB -n {scf_cores}
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "span[ptile={scf_ptile}]"

hostfile=`echo $LSB_DJOB_HOSTFILE`
NP=`cat $hostfile | wc -l`
cd $LS_SUBCWD
export PATH=/work/phy-huangj/app/vasp.5.4.4/bin:$PATH
COMMAND_std="mpirun -machinefile $LSB_DJOB_HOSTFILE -np $NP vasp_std"
'''

start_calc_command = 'bsub<'    #'sbatch' 'bsub<'