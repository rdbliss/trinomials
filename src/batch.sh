#!/bin/bash -l

# this is an example submission script for a slurm cluster. the exact
# parameters should be modified for your system.

#SBATCH --job-name=graphs

#SBATCH --array=2-500
#SBATCH -e ./errors/%a.txt
#SBATCH -o ./results/%a.txt

#SBATCH --cpus-per-task=40
#SBATCH --mem=2G

#SBATCH --time=1-00:00:00

module load python/3.13.2
module load flint
source /home/rwdb/graphs/venv/bin/activate
python ./graph.py $SLURM_ARRAY_TASK_ID
