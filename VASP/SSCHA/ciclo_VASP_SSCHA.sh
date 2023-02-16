#!/usr/bin/env bash
POPULATION=$1
SUPERCELL_SIZE=2
TEMPERATURE=300
NCONFSSCHA=512
NQIRR=4
echo "============================="
echo "Population="$POPULATION
echo "Supercell size="$SUPERCELL_SIZE
echo "Number of configurations="$NCONFSSCHA
echo "Temperature="$TEMPERATURE
echo "============================="
echo "Change directory to "pop$POPULATION
cd pop$POPULATION
echo "========"+`pwd`+"======="
cp ../POSCAR_UNITCELL POSCAR
python /home/diego/github/vasp-phonopy-sscha/vasp-phonopy-sscha/interface.py --f_processing $POPULATION $SUPERCELL_SIZE
python /home/diego/github/vasp-phonopy-sscha/vasp-phonopy-sscha/interface.py --en_processing $POPULATION
cd ..
##echo "Editing SSCHA inputfile"
##nano minimize.py
echo "------Running SSCHA-----"
#python3 minimize.py > minim$POPULATION.out
python3 minimnize.py -pop $POPULATION -nconf $NCONFSSCHA -cell $SUPERCELL_SIZE -temp $TEMPERATURE -nqirr $NQIRR > minim$POPULATION.out
echo "Change directory to "pop$(($POPULATION+1))
cd pop$(($POPULATION+1))
echo "========"+`pwd`+"======="
python /home/diego/github/vasp-phonopy-sscha/vasp-phonopy-sscha/interface.py --generate $NCONFSSCHA $(($POPULATION+1)) $SUPERCELL_SIZE $TEMPERATURE
cp ../POSCAR_UNITCELL POSCAR
python /home/diego/github/vasp-phonopy-sscha/vasp-phonopy-sscha/interface.py --to_vasp $(($POPULATION+1)) $SUPERCELL_SIZE
cp ../INCAR.sc vasp/INCAR
cp ../POTCAR.SrOTi vasp/POTCAR
cp ../ML_FF vasp/ML_FF
cp ../run.bash vasp/run.bash
nano vasp/run.bash
echo "---------------------------------------"
echo "Now is time to do the VASP calculations"
echo "---------------------------------------"
cd ..
