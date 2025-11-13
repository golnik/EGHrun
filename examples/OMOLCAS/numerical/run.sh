#!/bin/bash -l

#SBATCH --nodes 5
#SBATCH --ntasks-per-node 16
#SBATCH --cpus-per-task 1
#SBATCH --time=0:10:00
##SBATCH --nodelist=ltamp08

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module add gcc/8.5.0
module add openmpi/4.1.5/gcc8.5.0
module add anaconda/2024.02

source /software/apps/anaconda/2023.03/bin/activate
conda activate mpienv
export PATH="/software/apps/anaconda/2024.02/bin/:$PATH"
export PATH="/software/apps/openmpi/4.1.5/gcc8.5.0/bin/:$PATH"

geom_xyz="geometry.xyz"
geom_mol="geometry.mol"

script_template="molcas.in"

#tmp_dir="$PWD/calc"
tmp_dir="$TMPDIR/EGHrun/calc/"

run_out='output.out'
EGH_out='EGH.out'

EGHrun_path=~/programs/EGHrun/

time mpirun python $EGHrun_path/source/EGHrun.py -g $geom_xyz $geom_mol -rs $script_template -tdir $tmp_dir -run_out $run_out -EGH_out $EGH_out --print_script_out --calc_grad --calc_hess

exit 0
