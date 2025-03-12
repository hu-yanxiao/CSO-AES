# setup.py
from setuptools import setup, find_packages

setup(
    name='cso_aes',
    version='1.0',
    packages=find_packages(),
    author='Jing Huang',
    author_email='2760344463@qq.com',
    description='sus2_mlp_Activate_learning',
    install_requires=[
        'maml==2024.6.13',
        'matplotlib==3.8.0',
        'path==16.14.0',
        'ase==3.23.0',
        'dscribe==2.1.1',
        'Cython==3.0.11',
        'numpy==1.26.4',
    ]
)


'''
Remark:
install pymlip
install vaspkit
'ase==3.23.0',
'''
