units           metal
atom_style      atomic
boundary        p p p

read_data       tmp.lammps

pair_style      chimesFF
pair_coeff      * * params.txt

dump            1 all custom 1 XFSET id xu yu zu fx fy fz
dump_modify     1 format float "%20.15f"
run             0
