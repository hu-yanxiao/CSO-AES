from ase.io import iread
import sys

fin = sys.argv[1]
out = cfgs.replace('xyz', 'cfg')
map_dic={}
#ele = ['Li', 'B', 'O', 'Al', 'Si', 'P', 'S', 'Cl', 'Ga', 'Ge', 'As', 'Br', 'Sn', 'Sb', 'I']
#ele=['Ba', 'Pb', 'Ca', 'Sr', 'Bi', 'K', 'Na', 'Hf', 'Ti', 'Zr', 'Nb', 'Mg', 'In', 'Zn', 'O']
# ele = ['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'K', 'Ca',
#            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Rb',
#            'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',
#            'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
#            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Ac', 'Th', 'Pa', 'U', 'Np',
#            'Pu']
ele = ['O', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ca', 'Ti', 'V', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Zr', 'Nb', 'Mo', 'In', 'Sn', 'Sb', 'I', 'Hf', 'Ta', 'W', 'Bi']

for i in range(len(ele)):
    map_dic.update({ele[i]:i})
b=list(iread(fin))
ff = open(out,"w")
for i in range(len(b)):
    print(i)
    atoms=b[i]
    ele=atoms.get_chemical_symbols()
    nat=len(ele)
    cell = atoms.get_cell()
    pos = atoms.get_positions()
    force= atoms.get_forces()
    en=atoms.get_potential_energy()
    virial=atoms.info['virial']
    ff.write("""BEGIN_CFG\n""")
    ff.write(""" Size\n""")
    ff.write("""  {:6} \n""".format(nat))
    ff.write(""" Supercell \n""")
    ff.write("""{:15.10f} {:15.10f} {:15.10f}\n""".format(cell[0, 0], cell[0, 1], cell[0, 2]))
    ff.write("""{:15.10f} {:15.10f} {:15.10f}\n""".format(cell[1, 0], cell[1, 1], cell[1, 2]))
    ff.write("""{:15.10f} {:15.10f} {:15.10f}\n""".format(cell[2, 0], cell[2, 1], cell[2, 2]))
    ff.write("""AtomData:  id type       cartes_x      cartes_y      cartes_z     fx          fy          fz\n""")
    for i in range(nat):
        ff.write(
            """ {:6} {:6} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n""".format(i + 1, map_dic[ele[i]], pos[i, 0], pos[i, 1], pos[i, 2],
                                                                                         force[i,0], force[i,1], force[i,2]))
    ff.write("""Energy \n""")
    ff.write(f"""\t{en} \n""")
    ff.write("""PlusStress:  xx          yy          zz          yz          xz          xy \n""")
    ff.write(f"\t{virial[0,0]}  \t{virial[1,1]}  \t{virial[2,2]}  \t{virial[1,2]}  \t{virial[0,2]}  \t{virial[0,1]} \n")
    ff.write("""END_CFG \n""")
ff.close()
