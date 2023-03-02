#!/usr/bin/env bash
#
#  run_local.sh
#
#  Copyright 2023 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
POPULATION=$1       #population index
SUPERCELL_SIZE=2
TEMPERATURE=300
NCONFSSCHA=512      #number of configurations in the sscha ensemble
NQIRR=4             #number of irreducible q points
kong_liu_ratio=0.5  # The parameter that estimates whether the ensemble is still good
IONS=54             #number of atoms in the supercells
vasp_phonopy_sscha_interface_location="/media/diego/Calculations/SrTiO3/SSCHA/Hessian_vs_T_and_Nconfs/vasp-phonopy-sscha/vasp-phonopy-sscha"
echo "============================="
echo "Population="$POPULATION
echo "Supercell size="$SUPERCELL_SIZE
echo "Number of configurations="$NCONFSSCHA
echo "Temperature="$TEMPERATURE
echo "Number of atoms in the supercell="$IONS
echo "Number of irreducible q points="$NQIRR
echo "============================="
echo "Change directory to "pop$POPULATION
cd pop$POPULATION
echo "========"+`pwd`+"======="
cp ../POSCAR_UNITCELL POSCAR
python $vasp_phonopy_sscha_interface_location/interface.py --f_processing $POPULATION $SUPERCELL_SIZE
python $vasp_phonopy_sscha_interface_location/interface.py --en_processing $POPULATION
echo "============================="
echo "Change directory back"
cd ..
echo "========"+`pwd`+"======="
echo "------Running SSCHA-----"
python3 minimize.py -pop $POPULATION -nconf $NCONFSSCHA -cell $SUPERCELL_SIZE -temp $TEMPERATURE -nqirr $NQIRR > minim$POPULATION.out
echo "Change directory to "pop$(($POPULATION+1))
cd pop$(($POPULATION+1))
echo "========"+`pwd`+"======="
python $vasp_phonopy_sscha_interface_location/interface.py --generate $NCONFSSCHA $(($POPULATION+1)) $NQUIRR $TEMPERATURE
cp ../POSCAR_UNITCELL POSCAR
python $vasp_phonopy_sscha_interface_location/interface.py --to_vasp $(($POPULATION+1)) $SUPERCELL_SIZE
cp ../INCAR.sc vasp/INCAR
cp ../POTCAR.SrOTi vasp/POTCAR
cp ../ML_FF vasp/ML_FF
cp ../KPOINTS vasp/KPOINTS
mkdir vasp/forces
echo "============================="
echo "Change directory back"
cd ..
echo "========"+`pwd`+"======="
echo "----------------------------------------"
echo "Checking if the Kong-Liu parameter is OK"
echo "----------------------------------------"
kong_liu_1=`grep "Kong-Liu" minim1.out|head -1 | awk '{print $NF}'`
kong_liu_2=`grep "Kong-Liu" minim$POPULATION.out|tail -1 | awk '{print $NF}'`
echo "If this formula is OK, then you are converged:"
echo $(($kong_liu_1/$kong_liu_2))">"$kong_liu_ratio"?"
# if [[ $(($kong_liu_1/$kong_liu_2)) > $kong_liu_ratio ]]  #ALL THIS DOEN'T WORK!!
# then
#   echo "---------------------------------------"
#   echo "Now is time to do the VASP calculations"
#   echo "---------------------------------------"
#   echo "Change directory to "pop$(($POPULATION+1))"/vasp"
#   cd pop$(($POPULATION+1))/vasp
#   echo "========"+`pwd`+"======="
#   #########################################################################
#   #np=40           #number of cpus
#   #POPULATION=1    #population index
#   #IONS=54         #number of atoms in the supercells
#   #NCONFSSHA=300   #number of configurations in the sscha ensemble
#   #mpirun vasp_std > stdout
#   for i in `seq 1 $NCONFSSCHA`; do
#       echo -ne "---RUN--"$i"--of--"$NCONFSSCHA"\r"
#       echo `date` >> timing
#       cp POSCAR_$i POSCAR
#   #    mpirun -np $np vasp_std > stdout
#       ~/VASP/vasp.6.3.0/bin/vasp_std > stdout
#       grep "energy  without entropy" OUTCAR  >> energies
#       grep "forces" -A $IONS vasprun.xml > forces/forces_population$POPULATION'_'$i.dat
#       rm POSCAR
#       mv OUTCAR OUTCAR_$i
#       echo `date` >> timing
#   done
#   #########################################################################
# else
#   echo "--------------------------"
#   echo "The calculation converged."
#   echo "--------------------------"
# fi
echo "---------------------------------------"
echo "Now is time to do the VASP calculations"
echo "---------------------------------------"
echo "Change directory to "pop$(($POPULATION+1))"/vasp"
cd pop$(($POPULATION+1))/vasp
echo "========"+`pwd`+"======="
#########################################################################
for i in `seq 1 $NCONFSSCHA`; do
    echo -ne "---RUN--"$i"--of--"$NCONFSSCHA"\r"
    echo `date` >> timing
    cp POSCAR_$i POSCAR
    ~/VASP/vasp.6.3.0/bin/vasp_std > stdout
    grep "energy  without entropy" OUTCAR  >> energies
    grep "forces" -A $IONS vasprun.xml > forces/forces_population$(($POPULATION+1))'_'$i.dat
    rm POSCAR
    mv OUTCAR OUTCAR_$i
    echo `date` >> timing
done
