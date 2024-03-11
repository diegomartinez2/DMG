#!/bin/bash
#ORIGEN='/media/diego/Calculations/SrTiO3/SSCHA/Hessian_vs_T_and_Nconfs/VASP_17_4_2023/SrTiO3_3x3x3_50K_para_comparar_con_Espresso_Q4x4x4/'
ORIGEN='/media/diego/Calculations/SrTiO3/SSCHA/Hessian_vs_T_and_Nconfs/VASP_6_9_2023/4x4x4/16384/T300/'
#DESTINO='/scratch/diegom/SrTiO3/VASP/SrTiO3_3x3x3_50K_para_comparar_con_Espresso_Q4x4x4/'
DESTINO='/scratch/diegom/SrTiO3/VASP/SrTiO3_4x4x4_300K_16384'
CLUSTER='diegom@ekhi.cfm.ehu.es'
CHANGE="s/POPULATION=1/POPULATION=$1/g"
scp -r $ORIGEN"pop$1" $CLUSTER":"$DESTINO"/"
scp $ORIGEN"run_Ekhi.bash" $CLUSTER":"$DESTINO"/pop"$1"/vasp/run.sh"
#scp $ORIGEN"run_Ekhi1.bash" $CLUSTER":"$DESTINO"/pop"$1"/vasp/run1.sh"
#scp $ORIGEN"run_Ekhi2.bash" $CLUSTER":"$DESTINO"/pop"$1"/vasp/run2.sh"
#scp $ORIGEN"run_Ekhi3.bash" $CLUSTER":"$DESTINO"/pop"$1"/vasp/run3.sh"
#scp $ORIGEN"run_Ekhi4.bash" $CLUSTER":"$DESTINO"/pop"$1"/vasp/run4.sh"
#scp $ORIGEN"run_Ekhi5.bash" $CLUSTER":"$DESTINO"/pop"$1"/vasp/run5.sh"
#scp $ORIGEN"run_Ekhi6.bash" $CLUSTER":"$DESTINO"/pop"$1"/vasp/run6.sh"
#scp $ORIGEN"run_Ekhi7.bash" $CLUSTER":"$DESTINO"/pop"$1"/vasp/run7.sh"
#scp $ORIGEN"run_Ekhi8.bash" $CLUSTER":"$DESTINO"/pop"$1"/vasp/run8.sh"
#ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp;sed -i 's/POPULATION=1/POPULATION=$1/g' run.sh; sbatch run.sh"
echo "cd $DESTINO/pop$1/vasp; sbatch run.sh $1"
ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run.sh $1"
#ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run1.sh $1"
#sleep 60
#ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run2.sh $1"
#sleep 60
#ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run3.sh $1"
#sleep 60
#ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run4.sh $1"
#sleep 60
#ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run5.sh $1"
#sleep 60
#ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run6.sh $1"
#sleep 60
#ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run7.sh $1"
#sleep 60
#ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run8.sh $1"
#sleep 60


echo "##########"
while true
do
echo "Sleep..."
sleep 300
echo "Awake, checking:"
#files=$(ssh $CLUSTER 'squeue |grep diego')
#files=$(ssh $CLUSTER 'squeue -u diegom -n "16384_4_50_1,16384_4_50_2,16384_4_50_3,16384_4_50_4,16384_4_50_5,16384_4_50_6,16384_4_50_7,16384_4_50_8"|grep diego')
files=$(ssh $CLUSTER 'squeue -u diegom -n "16384_4_300"|grep diego')
if [[ $? != 0 ]]; then
    echo "No calculations running."
    echo "Hold for a while so we are sure the files are there"
    sleep 300
    echo "Taking the files from cluster"
    ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; rm ML_FF"
    scp -r $CLUSTER":"$DESTINO"/pop$1" $ORIGEN
    break
elif [[ $files ]]; then
    echo $files
else
    echo "No files found."
    #scp -r ${$CLUSTER":"$DESTINO"pop$1"} $ORIGEN
break
fi
done
echo "***DONE***"
