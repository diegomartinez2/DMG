for i in `seq 1 $N_RANDOM`
do
        cd ${i}
        sbash run.bash
done
