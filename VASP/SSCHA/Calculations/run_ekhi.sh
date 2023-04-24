#!/bin/bash
ORIGEN='/media/diego/Calculations/SrTiO3/SSCHA/Hessian_vs_T_and_Nconfs/VASP_17_4_2023/SrTiO3_3x3x3_50K_para_comparar_con_Espresso_Q4x4x4/'
DESTINO='/scratch/diegom/SrTiO3/VASP/SrTiO3_3x3x3_50K_para_comparar_con_Espresso_Q4x4x4/'
CLUSTER='diegom@ekhi.cfm.ehu.es'
scp -r $ORIGEN"pop$1" $CLUSTER$DESTINO
scp ${$ORIGEN"run_Ekhi.bash"} ${$CLUSTER":"$DESTINO"run.sh"}
ssh -t $CLUSTER "cd $DESTINO/pop$1/vasp; sbatch run.sh"
#scp -r /media/diego/Calculations/SrTiO3/SSCHA/Hessian_vs_T_and_Nconfs/VASP_17_4_2023/SrTiO3_3x3x3_50K_para_comparar_con_Espresso_Q4x4x4/pop$1 diegom@ekhi.cfm.ehu.es:/scratch/diegom/SrTiO3/VASP/SrTiO3_3x3x3_50K_para_comparar_con_Espresso_Q4x4x4/
#scp /media/diego/Calculations/SrTiO3/SSCHA/Hessian_vs_T_and_Nconfs/VASP_17_4_2023/SrTiO3_3x3x3_50K_para_comparar_con_Espresso_Q4x4x4/run_Ekhi.bash diegom@ekhi.cfm.ehu.es:/scratch/diegom/SrTiO3/VASP/SrTiO3_3x3x3_50K_para_comparar_con_Espresso_Q4x4x4/pop$1/vasp/run.sh
#ssh -t diegom@ekhi.cfm.ehu.es "cd /scratch/diegom/SrTiO3/VASP/SrTiO3_3x3x3_300K_para_comparar_con_Espresso_Q4x4x4/pop$1/vasp; sbatch run.sh"

##########
while True
do
echo "Sleep..."
sleep 300
echo "Awake, checking:"
files=$(ssh $CLUSTER 'squeue |grep diego')
#files=$(ssh diegom@ekhi.cfm.ehu.es 'squeue |grep diego')
if [[ $? != 0 ]]; then
    echo "Command failed."
elif [[ $files ]]; then
    echo $files
else
    echo "No files found."
    scp -r ${$CLUSTER":"$DESTINO"pop$1"} $ORIGEN
    #scp -r diegom@ekhi.cfm.ehu.es:/scratch/diegom/SrTiO3/VASP/SrTiO3_3x3x3_50K_para_comparar_con_Espresso_Q4x4x4/pop$1 /media/diego/Calculations/SrTiO3/SSCHA/Hessian_vs_T_and_Nconfs/VASP_17_4_2023/SrTiO3_3x3x3_50K_para_comparar_con_Espresso_Q4x4x4/
break
fi
done
