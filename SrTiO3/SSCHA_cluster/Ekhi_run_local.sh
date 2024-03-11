#!/bin/bash
#ORIGEN='/media/diego/Calculations/SrTiO3/SSCHA/Hessian_vs_T_and_Nconfs/VASP_6_9_2023/4x4x4/16384/T300/'
#DESTINO='/scratch/diegom/SrTiO3/VASP/SrTiO3_4x4x4_300K_16384'
#CLUSTER='diegom@ekhi.cfm.ehu.es'
#CHANGE="s/POPULATION=1/POPULATION=$1/g"
#scp -r $ORIGEN"pop$1" $CLUSTER":"$DESTINO"/"
#scp $ORIGEN"run_Ekhi.bash" $CLUSTER":"$DESTINO"/pop"$1"/vasp/run.sh"
#echo "cd $DESTINO/pop$1/vasp; sbatch run.sh $1"
#ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run.sh $1"
sbatch run_Ekhi.bash

echo "##########"
while true
do
echo "Sleep..."
sleep 300
echo "Awake, checking:"
#files=$(ssh $CLUSTER 'squeue |grep diego')
#files=$(ssh $CLUSTER 'squeue -u diegom -n "16384_4_50_1,16384_4_50_2,16384_4_50_3,16384_4_50_4,16384_4_50_5,16384_4_50_6,16384_4_50_7,16384_4_50_8"|grep diego')
files=$(ssh $CLUSTER 'squeue -u diegom -n "tetragonal_2"|grep diego')
if [[ $? != 0 ]]; then
    echo "No calculations running."
    echo "Hold for a while so we are sure the files are there"
    sleep 300
    #echo "Taking the files from cluster"
    #ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; rm ML_FF"
    #scp -r $CLUSTER":"$DESTINO"/pop$1" $ORIGEN
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
