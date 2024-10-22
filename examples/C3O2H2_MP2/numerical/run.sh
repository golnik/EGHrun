#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --ntasks-per-node 10
#SBATCH --cpus-per-task 1
#SBATCH --time=0:10:00
#SBATCH --account=lcpt
#SBATCH --exclusive
#SBATCH --mem 8000mb

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load intel/18.0.5
module load intel-mpi/2018.4.274
module load python/3.7.3

source /home/ngolubev/Packages/virtualenv/x86_E5v4_Mellanox_intel/bin/activate

geom_xyz="geometry.xyz"
geom_mol="geometry.mol"

script_template="scf_adc.in"

tmp_dir="$TMPDIR/EGHrun/calc/"

run_out='output.out'
EGH_out='EGH.out'

time srun python /home/ngolubev/programs/EGHrun/source/EGHrun.py -g $geom_xyz $geom_mol -rs $script_template -tdir $tmp_dir --calc_grad --calc_hess --z_sym -run_out $run_out -EGH_out $EGH_out

exit 0
