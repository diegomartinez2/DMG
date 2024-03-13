#!/bin/bash
echo "cd pop$1/vasp; sbatch run.sh $1"
cd "pop$1/vasp";pwd
cp ../../run_Ekhi.bash run_Ekhi.bash
sbatch run_Ekhi.bash
echo "return to *base* directoty"
cd ../..; pwd
echo "##########"
while true
do
echo "Sleep..."
sleep 300
echo "Awake, checking:"
#files=$(ssh $CLUSTER 'squeue -u diegom -n "16384_4_50_1,16384_4_50_2,16384_4_50_3,16384_4_50_4,16384_4_50_5,16384_4_50_6,16384_4_50_7,16384_4_50_8"|grep diego')
files=$(squeue -u diegom -n "tetragonal_2"|grep diego)
if [[ $? != 0 ]]; then
    echo "No calculations running."
    echo "Hold for a while so we are sure the files are there"
    sleep 300
    break
elif [[ $files ]]; then
    echo $files
else
    echo "No files found."
break
fi
done
echo "***DONE***"
