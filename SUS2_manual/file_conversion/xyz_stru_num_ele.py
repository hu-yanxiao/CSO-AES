import sys
from ase.io import iread,read
from ase.data import atomic_numbers

input = sys.argv[1]
data = iread(input)

count = 0
temp = set()
for a in data:
    count += 1
    temp = temp | set(a.get_chemical_symbols())
print(f'stru_num: {count}')
print(sorted(temp, key=lambda x: atomic_numbers[x]))
print(f'num_ele: {len(temp)}')