#! /bin/bash 
#SBATCH --job-name=temp
#SBATCH --mem-per-cpu=2000mb
#SBATCH --time=01:00:00
#SBATCH --output=prod_%A-%a.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1            # Number of cores per MPI rank 
#SBATCH --ntasks-per-node=1         # How many tasks on each node
#SBATCH --ntasks-per-socket=1        # How many tasks on each CPU or socket
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically on nodes and sockets
#SBATCH --array=0-1
#SBATCH  --account=alberto.perezant
#SBATCH  --qos=alberto.perezant




b=`perl -e 'printf("%02i", $ARGV[0]);' ${SLURM_ARRAY_TASK_ID}`
extract_trajectory extract_traj_dcd --replica ${SLURM_ARRAY_TASK_ID} trajectory.$b.dcd
#extract_trajectory follow_dcd --replica ${SLURM_ARRAY_TASK_ID} follow.$b.dcd
