#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --time=1:00:00
#SBATCH --account=theoryexpt
#SBATCH --partition=standard

module load anaconda
source /opt/ohpc/pub/apps/anaconda/2024.06/bin/activate
conda activate oceloteMpiEnv
python PES.py

exit 0
