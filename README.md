# CSO-AES
Covering Set Optimazation driven Atomic Environment Sampling strategy for machine learning interatomic potential   

# Installation
conda create --name cso python=3.10  
cd CSO-AES/cso_aes  
python setup.py install  
pip install -r requirement.txt  
if scf_cal_engine = ABACUS, please  install ase-abacus (pip install git https://gitlab.com/1041176461/ase-abacus.git).  
if scf_cal_engine = VASP, please install [VASPKIT](https://vaspkit.com/installation.html) to automatically generate POTCAR file.

# Additional Installation Dependencies for the Program
# 1. SUS2-MLIP installation  
install [SUS2-MLIP](https://github.com/hu-yanxiao/SUS2-MLIP)
# 2. pysus2mlip installation
tar -zxvf pysus2mlip.tar.gz  
cd pysus2mlip  
Modify the correct SUS2-related file locations as shown in the content of setup.py. Then execute: "CC=icc CXX=icpc pip install -e ."    

![pysus2mlip](https://github.com/user-attachments/assets/c432dc3b-16b2-4ac0-8e53-f3ade9aec096)

# 3. LAMMPS installation
SUS2-MLIP models can be used in [LAMMPS](https://github.com/lammps/lammps) simulation via the interface [interface-lammps-mlip-v2](https://gitlab.com/ashapeev/interface-lammps-mlip-2/-/tree/master?ref_type=heads).

# Tutorial
Video: [Install](https://b23.tv/PSIvqp5) [Run](https://b23.tv/dy2E1WQ) [introduction](https://b23.tv/2xEr8Tz)
Sample Manual: [Onenote](https://pan.quark.cn/s/85215f61ed71?pwd=r8Fz)
