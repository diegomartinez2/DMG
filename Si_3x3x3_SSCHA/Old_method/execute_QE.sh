N_RANDOM=`grep N_RANDOM create_configurations.py | head -1 | awk '{print $3}'`
POPULATION=`grep POPULATION create_configurations.py | head -1 | awk '{print $3}'`

# Extract forces, energies, and stresses

cp harmonic_calculation/*upf population${POPULATION}_ensemble
cd population${POPULATION}_ensemble

for i in `seq 1 $N_RANDOM`
do
        mpirun -np 8 pw.x < scf_${i}.in > scf_${i}.out
        rm -rf tmp
done

