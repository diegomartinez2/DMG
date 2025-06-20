#!/bin/bash
echo "Start this with phonopy"
echo "this creates the SPOSCAR and POSCAR-* files for a 2x2x2 supercell, make it bigger for better convergence of the DOS"
phonopy -c POSCAR --dim="2 2 2" -d
echo "...continue with atomsk"
atomsk POSCAR-001 lammps
atomsk POSCAR-002 lammps
atomsk POSCAR-003 lammps
atomsk POSCAR-004 lammps
atomsk POSCAR-005 lammps
atomsk POSCAR-006 lammps
atomsk POSCAR-007 lammps
atomsk POSCAR-008 lammps
atomsk POSCAR-009 lammps
atomsk POSCAR-010 lammps
atomsk POSCAR-011 lammps
atomsk POSCAR-012 lammps
atomsk POSCAR-013 lammps
atomsk POSCAR-014 lammps
atomsk POSCAR-015 lammps
atomsk POSCAR-016 lammps
atomsk POSCAR-017 lammps
atomsk POSCAR-018 lammps
atomsk POSCAR-019 lammps
atomsk POSCAR-020 lammps
atomsk POSCAR-021 lammps
atomsk POSCAR-022 lammps
atomsk POSCAR-023 lammps
atomsk POSCAR-024 lammps
atomsk POSCAR-025 lammps
atomsk POSCAR-026 lammps
atomsk POSCAR-027 lammps
atomsk POSCAR-028 lammps
atomsk POSCAR-029 lammps
atomsk POSCAR-030 lammps
atomsk POSCAR-031 lammps
atomsk POSCAR-032 lammps
atomsk POSCAR-033 lammps
atomsk POSCAR-034 lammps
atomsk POSCAR-035 lammps
atomsk POSCAR-036 lammps
atomsk POSCAR-037 lammps
atomsk POSCAR-038 lammps
atomsk POSCAR-039 lammps
atomsk POSCAR-040 lammps
atomsk POSCAR-041 lammps
atomsk POSCAR-042 lammps
atomsk POSCAR-043 lammps
atomsk POSCAR-044 lammps
atomsk POSCAR-045 lammps
#for i in POSCAR-*; do
#  atomsk $i lammps
#done  
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
