#!/bin/bash
[ -n "$DEBUG" ] && set -x

TASK[0]=0
export TASK

task1() {
    file="/tmp/task1_progress"; echo 0 > $file
    for i in $(seq 0 10 100); do echo $i > $file; sleep 1; done
}

task2() {
    file="/tmp/task2_progress"; echo 0 > $file
    for i in $(seq 0 10 100); do echo $i > $file; sleep 2; done
}

task1 & task2 &

unset TIME_TO_FINISH PROGRESS STATUS
while /bin/true;
do
    unset TASKS
    for i in 1 2;
    do
        PROGRESS[$i]=$(cat /tmp/task${i}_progress)
        if [ ${PROGRESS[$i]} -eq 100 ];
        then
            STATUS[$i]=0
        else
            STATUS[$i]=-${PROGRESS[$i]}
        fi

        TASKS+=("Task $i" "${STATUS[$i]}")
    done

    # 0: success
    # 1: failed
    # 2: passed
    # 3: completed
    # 4: checked
    # 5: done
    # 6: skipped
    # 7: in progress
    # -X: 0-100, progress of process
    dialog \
        --title "Mixed gauge demonstration" \
        --backtitle "Backtitle" \
        --mixedgauge "This is a prompt message,\nand this is the second line." \
            0 0 $(($((${PROGRESS[1]}+${PROGRESS[2]}))/2)) \
            "${TASKS[@]}"

    I=0
    for a in 1 2; do [ ${PROGRESS[$a]} -ge 100 ] && I=$((I+1)); done

    [ $I -eq 2 ] && break
    sleep 1
done
echo waiting
wait

set +x
