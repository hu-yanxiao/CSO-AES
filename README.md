# CSO-AES
Covering Set Optimazation driven Atomic Environment Sampling strategy for machine learning interatomic potential 

# Installationpysus2mlip.tar.gz
conda create --name cso python=3.10

cd CSO-AES/cso_aes

python setup.py install

pip install -r requirement.txt

if scf_cal_engine = ABACUS, please  install ase-abacus (pip install git https://gitlab.com/1041176461/ase-abacus.git) .

if scf_cal_engine = VASP, please install VASPKIT (https://vaspkit.com/installation.html) to automatically generate POTCAR file.


# pysus2mlip installation
tar -zxvf pysus2mlip.tar.gz
cd pysus2mlip

modify the location of the sus2-related files in the program. Then execute: CC=icc CXX=icpc pip install -e . 

![pysus2mlip](https://github.com/user-attachments/assets/c432dc3b-16b2-4ac0-8e53-f3ade9aec096)

