import os
import numpy as np
from ase.io import read,write,iread


# Constants definition
HA_TO_EV = 27.211386245988  # 1 Ha = 27.211 eV
HA_BOHR_TO_EV_ANGSTROM = 51.422086190832  # 1 Ha/Bohr = 51.422 eV/Å
GPA_TO_EV_A3 = 0.0062415091258833  # 1 GPa = 0.00624 eV/Å³

def parse_forces_and_stress(logout_path):
    """ Parse the logout file to extract energy, atomic forces, and stress tensor """
    energy, forces, stress = None, [], None
    with open(logout_path, 'r') as f:
        lines = f.readlines()

    # Extract energy
    for line in lines:
        if 'ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]' in line:
            energy = float(line.split()[-1]) * HA_TO_EV  # Convert to eV
            break

    # Extract atomic forces
    force_start, force_end = None, None
    for i, line in enumerate(lines):
        if 'ATOMIC FORCES in [a.u.]' in line:
            force_start = i + 2
        if 'SUM OF ATOMIC FORCES' in line:
            force_end = i
            break
    forces = [
        [float(x) * HA_BOHR_TO_EV_ANGSTROM for x in line.split()[3:6]]
        for line in lines[force_start:force_end]
        if line.split()[0].isdigit()
    ]

    # Extract stress tensor
    for i, line in enumerate(lines):
        if 'STRESS| Analytical stress tensor [GPa]' in line:
            stress_lines = lines[i + 2:i + 5]
            stress = np.array([
                [float(x) for x in row.split()[2:]]
                for row in stress_lines
            ]) * GPA_TO_EV_A3  # Convert to eV/Å³
            break

    return energy, forces, stress


def calculate_virial(stress, lattice):
    """ Calculate the virial and return it as a string """
    volume = np.linalg.det(lattice)  # Calculate unit cell volume
    virial = (stress.flatten() * volume)  # Stress tensor multiplied by volume
    return ' '.join(f'{v:.6f}' for v in virial)


def update_cp2k_file(cp2k_xyz_path, energy, forces, stress, new_file_path):
    """ Update the cp2k.xyz file with energy, forces, and stress tensor """
    with open(cp2k_xyz_path, 'r') as f:
        lines = f.readlines()

    # Extract lattice constants and compute volume
    lattice_str = lines[1].split('Lattice=')[1].split('Properties')[0].strip().replace('"', '')
    lattice = np.array([list(map(float, lattice_str.split()))[i:i + 3] for i in range(0, 9, 3)])

    # Calculate virial and replace the virial entry
    virial_str = calculate_virial(stress, lattice)

    # Update the second line with energy, forces, and stress tensor
    for i, line in enumerate(lines):
        if i == 1:
            energy_str = line.split('energy=')[1].split()[0]
            new_line = line.replace(f'energy={energy_str}', f'energy={energy:.6f}')
            start_idx = new_line.find('virial="') + len('virial="')
            end_idx = new_line.find('"', start_idx)
            new_line = new_line[:start_idx] + virial_str + new_line[end_idx:]
            lines[i] = new_line
        elif i >= 2 and i < 2 + len(forces):  # Update atomic forces from line 3 onward
            parts = line.split()
            # Format each column to align properly: element, position (x, y, z), force (fx, fy, fz)
            formatted_line = f"{parts[0]:<2}  " + \
                             f"{float(parts[1]):>15.8f}  {float(parts[2]):>15.8f}  {float(parts[3]):>15.8f}  " + \
                             f"{forces[i - 2][0]:>15.8f}  {forces[i - 2][1]:>15.8f}  {forces[i - 2][2]:>15.8f}\n"
            lines[i] = formatted_line

    # Write the updated content to the new file
    with open(new_file_path, 'w') as f:
        f.writelines(lines)

def collect_efs(input_path):
    xyz = os.path.join(input_path,'cp2k.xyz')
    logout_path = os.path.join(input_path, 'logout')
    energy, forces, stress = parse_forces_and_stress(logout_path)
    old_atoms = read(xyz)
    temp_poscar = os.path.join(input_path,'temp.vasp')
    write(temp_poscar,old_atoms ,format='vasp')
    atoms = read(temp_poscar)
    atoms.info['energy'] = energy
    atoms.info['virial'] = np.array(stress)
    atoms.info['pbc'] = "T T T"
    atoms.arrays['forces'] = np.array(forces)
    atoms.pbc = [True, True, True]
    if os.path.exists(temp_poscar):
        os.remove(temp_poscar)
    return atoms


if __name__ == '__main__':
    input_path = os.getcwd()
    atoms = cp2k2xyz(input_path)
    print(atoms)
    write('test.xyz',atoms,format='extxyz')

    # """ Main function to read logout, process data, and update cp2k.xyz """
    # energy, forces, stress = parse_forces_and_stress(logout_path)
    # update_cp2k_file(cp2k_xyz_path, energy, forces, stress, 'cp2knew.xyz')
