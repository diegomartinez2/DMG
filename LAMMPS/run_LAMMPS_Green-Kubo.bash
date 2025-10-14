for seed in 12345 67890 45678; do
    sed "s/12345/${seed}/g" kappa.in > kappa_${seed}.in
    lmp_serial -in kappa_${seed}.in -log log_${seed}.lammps
done
