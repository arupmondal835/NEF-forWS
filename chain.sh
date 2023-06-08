#!/bin/bash

script=MELD_job.sh
njob=4

echo I\'ll submit the script $script will restart $njob times
echo If this is not what you want you have 10s to kill this process
echo The positional arguments are: PBS_script name number_of_restarts
#sleep 10


PID=$(sbatch $script|awk '{print $NF}')
echo $PID
for i in `seq 1 $njob`; do
        PID=$(sbatch --dependency=afterany:$PID  $script| awk '{print $NF}')
        echo $PID
    done

script=/orange/alberto.perezant/arup.mondal/NEF/test_for_web_server/array_job_temp.sh


PID_3=$(sbatch --dependency=afterany:$PID  $script| awk '{print $NF}')
echo $PID_3

script=/orange/alberto.perezant/arup.mondal/NEF/test_for_web_server/hierarchical_clustering.sh 

PID_4=$(sbatch --dependency=afterany:$PID_3  $script| awk '{print $NF}')
echo $PID_4


