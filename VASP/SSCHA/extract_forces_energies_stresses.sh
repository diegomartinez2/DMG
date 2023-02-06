N_RANDOM=`grep N_RANDOM create_configurations.py | head -1 | awk '{print $3}'`
POPULATION=`grep POPULATION create_configurations.py | head -1 | awk '{print $3}'`
#NAT=`grep 'number of atoms/cell      =' 'population'$POPULATION'_ensemble/scf_1.out' | awk '{print $5}'`
NAT=`grep 'NIONS' 'population'$POPULATION'_ensemble/1/OUTCAR' | awk '{print $12}'`
# Define constants
Ev_to_Ry=0.073498644351
Angst_to_bohr=1.8897259886
# Extract forces, energies, and stresses (CHANGE IT TO WORK WITH THE VASP OUTPUT)
#Notes: ,1 eV = 0.073498644351 Ry

for i in `seq 1 $N_RANDOM`
do
#        grep force 'population'$POPULATION'_ensemble/scf_'$i'.out' | grep atom | head -${NAT} | awk '{print $7,$8,$9}' > 'population'$POPULATION'_ensemble/forces_population'$POPULATION'_'$i'.dat'
#        grep -A3 'total   stress' 'population'$POPULATION'_ensemble/'$i'/OUTCAR' | tail -3 | awk '{print $1,$2,$3}' > 'population'$POPULATION'_ensemble/pressures_population'$POPULATION'_'$i'.dat'
        grep -A3 '=-STRESS' 'population'$POPULATION'_ensemble/'$i'/OUTCAR' | tail -1 | awk '{print $2,$5,$7"\n"$5,$2,$6"\n"$7,$6,$3}' > 'population'$POPULATION'_ensemble/pressures_population'$POPULATION'_'$i'.dat'
done

#for i in `seq 1 $N_RANDOM`; do grep ! 'population'$POPULATION'_ensemble/scf_'$i'.out' | awk '{print $5}'; done > 'population'$POPULATION'_ensemble/energies_supercell_population'$POPULATION'.dat'
for i in `seq 1 $N_RANDOM`
do
  grep "free  energy ML TOTEN  =" 'population'$POPULATION'_ensemble/'$i'/OUTCAR' | tail -1 | awk '{printf "%6.6f \n ", $6 * 0.0734986176 }'
done > 'population'$POPULATION'_ensemble/energies_supercell_population'$POPULATION'.dat'

for i in `seq 1 $N_RANDOM`
do
  N_ions=$(grep NIONS 'population'$POPULATION'_ensemble/'$i'/OUTCAR' |awk '{printf"%d",$NF}')
  awk '/TOTAL-FORCE \(eV\/Angst\)/{
	n++
	{
		printf("\nN_iteration : %5d\n",n)
		printf("%5s %10s %10s %10s %12s %12s %12s\n","atom","X","Y","Z","Fx","Fy","Fz")
		for(i=1;i<='${N_ions}';i++){
			getline
			{
				printf("%5d %10.5f %10.5f %10.5f %12.6f %12.6f %12.6f\n",i,$1*Angst_to_bohr,$2*Angst_to_bohr,$3*Angst_to_bohr,$4*(Ev_to_Ry/Angst_to_bohr),$5*(Ev_to_Ry/Angst_to_bohr),$6*(Ev_to_Ry/Angst_to_bohr))
			}
		}

	}
}' 'population'$POPULATION'_ensemble/'$i'/OUTCAR' |tee  'population'$POPULATION'_ensemble/forces_population'$POPULATION'_'$i'.dat'
