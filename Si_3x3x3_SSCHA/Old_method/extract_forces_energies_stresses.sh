N_RANDOM=`grep N_RANDOM create_configurations.py | head -1 | awk '{print $3}'`
POPULATION=`grep POPULATION create_configurations.py | head -1 | awk '{print $3}'`
NAT=`grep 'number of atoms/cell      =' 'population'$POPULATION'_ensemble/scf_1.out' | awk '{print $5}'`

# Extract forces, energies, and stresses

for i in `seq 1 $N_RANDOM`
do
        grep force 'population'$POPULATION'_ensemble/scf_'$i'.out' | grep atom | head -${NAT} | awk '{print $7,$8,$9}' > 'population'$POPULATION'_ensemble/forces_population'$POPULATION'_'$i'.dat'
        grep -A3 'total   stress' 'population'$POPULATION'_ensemble/scf_'$i'.out' | tail -3 | awk '{print $1,$2,$3}' > 'population'$POPULATION'_ensemble/pressures_population'$POPULATION'_'$i'.dat'
done

for i in `seq 1 $N_RANDOM`; do grep ! 'population'$POPULATION'_ensemble/scf_'$i'.out' | awk '{print $5}'; done > 'population'$POPULATION'_ensemble/energies_supercell_population'$POPULATION'.dat'

