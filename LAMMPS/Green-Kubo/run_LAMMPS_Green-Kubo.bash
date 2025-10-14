for seed in 12345 23456 34567; do
    sed "s/12345/${seed}/g" kappa.in > kappa_${seed}.in
    lmp_serial -in kappa_${seed}.in -log log_${seed}.lammps
done
