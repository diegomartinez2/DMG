POPULATION=1
SUPERCELL_SIZE=4
TEMPERATURE=300
NCONFSSCHA=512

cp ../POSCAR_UNITCELL POSCAR
python /home/diego/github/vasp-phonopy-sscha/vasp-phonopy-sscha/interface.py --f_processing $POPULATION $SUPERCELL_SIZE
python /home/diego/github/vasp-phonopy-sscha/vasp-phonopy-sscha/interface.py --en_processing $POPULATION
cd ..
nano minimize.py
python3 minimize.py > minim$POPULATION.out
cd pop{&POPULATION+1}
python /home/diego/github/vasp-phonopy-sscha/vasp-phonopy-sscha/interface.py --generate $NCONFSSCHA {&POPULATION+1} $SUPERCELL_SIZE $TEMPERATURE
