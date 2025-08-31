import numpy as np
import sys
map_dic={}
#ele = ['Li', 'B', 'O', 'Al', 'Si', 'P', 'S', 'Cl', 'Ga', 'Ge', 'As', 'Br', 'Sn', 'Sb', 'I']
#ele = ['H', 'O', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Sr', 'La']
# ele = ['Li', 'Be', 'Na', 'Mg', 'Al', 'Si', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
#        'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
#        'In', 'Sn', 'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Lu', 'Hf',
#        'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Pb']

ele = ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ca',
       'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Rb',
       'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',
       'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
       'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Ac', 'Th', 'Pa',
       'U', 'Np', 'Pu']

for i in range(len(ele)):
    map_dic.update({i:ele[i]})
print(map_dic)


cfgs = sys.argv[1]
out = cfgs.replace('cfg', 'xyz')
# cfgs="./train.cfg"
# out = 'NaCaTiO3.xyz'

with open(cfgs) as f:
    lines = f.readlines()
    cfgcnt = 0
    for line in lines:
        if line == ' Size\n':
            cfgcnt += 1

    cntr=1
    for i in range(len(lines)):
        if lines[i] != 'BEGIN_CFG\n':
            continue
#        else:
#            print("reading cfg#"+str(cntr+1))
        size = int(lines[i+2].split()[0])
        energy=float(lines[i+9+size].split()[0])
        lat=lines[i+4].split()+lines[i+5].split()+lines[i+6].split()
        for k in range(len(lat)):
            lat[k]=float(lat[k])
        stress=lines[i+11+size].split()
        for l in range(len(stress)):
            stress[l]=float(stress[l])
        _stress=[stress[0],stress[5],stress[4],stress[5],stress[1],stress[3],stress[4],stress[3],stress[2]]


        nruter = []
        for j in range(size):
            tmp = []
            words = lines[i+8+j].split()
            tmp.append(map_dic[int(words[1])])
            tmp.append(float(words[2]))
            tmp.append(float(words[3]))
            tmp.append(float(words[4]))
            tmp.append(float(words[-3]))
            tmp.append(float(words[-2]))
            tmp.append(float(words[-1]))
            nruter.append(tmp)


        with open(out,'a') as ff:
            ff.write(str(size)+'\n')
            ff.write("""Lattice="{: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} " Properties=species:S:1:pos:R:3:forces:R:3 energy={: 12.8f} virial= "{: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f}" pbc="T T T" \n""".format(*lat,energy,*_stress))
            for k in range (size):
                ff.write('{} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f} {: 12.8f}\n'.format(*nruter[k]))
        print(cntr)
        cntr += 1

