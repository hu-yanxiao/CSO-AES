from ase.io import iread,read,write
from ase.data import chemical_symbols

def xyz2cfg(f_in):
    atoms_list = list(iread(f_in))
    ele_set = set()
    ele_list = []
    for a in atoms_list:
        temp_list = a.get_chemical_symbols()
        ele_list += temp_list
        temp_set = set(temp_list)
        ele_set = temp_set | ele_set
    elements = list(ele_set)
    sorted_elements = sorted(elements, key=lambda x: chemical_symbols.index(x))

    map_dic = {element: i for i, element in enumerate(sorted_elements)}
    print(map_dic)
    print(list(map_dic.keys()))

    b = list(iread(f_in))
    ff = open(f_in.replace('xyz', 'cfg'), "w")

    for i in range(len(b)):
        #print(i)
        atoms = b[i]
        ele = atoms.get_chemical_symbols()
        nat = len(ele)
        cell = atoms.get_cell()
        pos = atoms.get_positions()
        force = atoms.get_forces()
        en = atoms.get_potential_energy()
        #label = atoms.info['label']
        virial = atoms.info['virial']
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
                """ {:6} {:6} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n""".format(i + 1, map_dic[ele[i]],
                                                                                                pos[i, 0], pos[i, 1],
                                                                                                pos[i, 2],
                                                                                                force[i, 0],
                                                                                                force[i, 1],
                                                                                                force[i, 2]))
        ff.write("""Energy \n""")
        ff.write(f"""\t{en} \n""")
        ff.write("""PlusStress:  xx          yy          zz          yz          xz          xy \n""")
        ff.write(
            f"\t{virial[0, 0]}  \t{virial[1, 1]}  \t{virial[2, 2]}  \t{virial[1, 2]}  \t{virial[0, 2]}  \t{virial[0, 1]} \n")
        #ff.write(f"""Feature   uid {label} \n""")
        ff.write("""END_CFG \n""")
    ff.close()

if __name__ =='__main__':
    import sys
    f_in = sys.argv[1]
    xyz2cfg(f_in)