#!/bin/bash
echo "Start this with phonopy"
echo "this creates the SPOSCAR and POSCAR-* files for a 2x2x2 supercell, make it bigger for better convergence of the DOS"
phonopy -c POSCAR --dim="2 2 2" -d
echo "...continue with atomsk"
for i in POSCAR-*; do
  atomsk $i lammps
done
echo "...now convert to LAMMPS with charges"
python Convert_to_LAMMPS_charges.py
echo "...calculate with LAMMPS"
for i in POSCAR-*.lammps; do
    base=$(basename $i .lammps)
    /home/rdettori/lmp_mpi_chimes -in in.phonopy -v datafile $i -v dumpfile forces-${base}.dump
done
echo "...Taking the forces to phonopy"
phonopy -f forces-POSCAR-*.dump
echo "...calculate DOS. The mesh of 20x20x20 is reasonable, but increase to 30x30x30 for smother results"
phonopy -c POSCAR --dim="2 2 2" --mesh="20 20 20" -t
phonopy -p mesh.yaml -s
echo "...alt..."
phonopy -c POSCAR -f forces-POSCAR-*.dump --dim="2 2 2" --band="band.conf"
phonopy -p band.yaml -s

echo "---END---"
