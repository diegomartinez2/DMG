N_RANDOM=`grep N_RANDOM create_configurations.py | head -1 | awk '{print $3}'`

for i in `seq 1 $N_RANDOM`
do
#        ln -sf ./ML_FF ./${i}/ML_FF
        cd ${i}
        sbash run.bash
done
