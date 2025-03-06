
npt = '''variable dt equal 0.001
variable StepNPT equal 100000
variable PreStep equal 0
variable nevery equal 100
variable out_dump_file_0 string force.0.dump
# variable out_dump_file string
# variable mlip_ini string

units metal
boundary p p p
atom_style atomic

box tilt large
read_data data.in


pair_style  mlip ${mlip_ini}
pair_coeff * *

neighbor 2.0 bin
neigh_modify delay 10 check yes

timestep    ${dt}
variable Tdamp equal "v_dt * 100"
variable Pdamp equal "v_dt * 1000"
velocity        all create ${T} ${random} mom yes rot yes dist gaussian
thermo 100

fix 100 all npt temp ${T} ${T} ${Tdamp} aniso 0.0 0.0 ${Pdamp}
velocity all scale ${T}
velocity all zero linear

run             ${PreStep}

reset_timestep 0

compute pe all pe/atom
dump 0 all custom ${nevery} ${out_dump_file_0} id type x y z fx fy fz c_pe
dump_modify 0 format line "%6d %3d %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e" sort id

dump 1 all custom ${nevery} ${out_dump_file} id type x y z fx fy fz c_pe
dump_modify 1 sort id

run             ${StepNPT}
'''

nvt = '''variable dt equal 0.001
variable StepNVT equal 100000
variable PreStep equal 0
variable nevery equal 100
variable out_dump_file_0 string force.0.dump
# variable out_dump_file string
# variable mlip_ini string

units metal
boundary p p p
atom_style atomic

box tilt large
read_data data.in


pair_style  mlip ${mlip_ini}
pair_coeff * *

neighbor 2.0 bin
neigh_modify delay 10 check yes

timestep    ${dt}
variable Tdamp equal "v_dt * 100"
variable Pdamp equal "v_dt * 1000"
velocity        all create ${T} ${random} mom yes rot yes dist gaussian
thermo 100

fix 100 all nvt temp ${T} ${T} ${Tdamp}
velocity all scale ${T}
velocity all zero linear

run             ${PreStep}

reset_timestep 0

compute pe all pe/atom
dump 0 all custom ${nevery} ${out_dump_file_0} id type x y z fx fy fz c_pe
dump_modify 0 format line "%6d %3d %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e" sort id

dump 1 all custom ${nevery} ${out_dump_file} id type x y z fx fy fz c_pe
dump_modify 1 sort id

run             ${StepNVT}
'''