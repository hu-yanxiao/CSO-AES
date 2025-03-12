# CSO-AES
Covering Set Optimazation driven Atomic Environment Sampling strategy for machine learning interatomic potential 

# Installation
conda create --name cso python=3.10

cd CSO-AES/cso_aes

python setup.py install

if scf_cal_engine = ABACUS, please  install ase-abacus (pip install git https://gitlab.com/1041176461/ase-abacus.git) .

if scf_cal_engine = VASP, please install VASPKIT (https://vaspkit.com/installation.html) to automatically generate POTCAR file.
