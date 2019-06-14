#!/bin/bash -l

#SBATCH --nodes 2
#SBATCH --ntasks-per-node 7
#SBATCH --cpus-per-task 4
#SBATCH --time=1:00:00
#SBATCH --account=lcpt
#SBATCH --exclusive
#SBATCH --mem 8000mb

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load intel/18.0.2
module load intel-mpi/2018.2.199
module load python/3.6.5

geom_xyz="geometry.xyz"
geom_mol="geometry.mol"

script_template="scf_adc.in"

tmp_dir="$TMPDIR/EGHrun/calc/"

run_out='output.out'
EGH_out='EGH.out'

#activate python environment and run calculations
source /home/ngolubev/Packages/virtualenv/x86_E5v4_Mellanox_intel/bin/activate
time srun python ../../source/EGHrun.py -g $geom_xyz $geom_mol -rs $script_template -tdir $tmp_dir --calc_grad --z_sym

exit 0
