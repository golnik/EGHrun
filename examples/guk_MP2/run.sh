#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --ntasks-per-node 11
#SBATCH --cpus-per-task 1
#SBATCH --time=0:10:00
#SBATCH --nodelist=ltamp08

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

script_template="scf_mp2.in"

#tmp_dir=$PWD/calc
#tmp_dir="$SLURM_SUBMIT_DIR/EGHrun/calc/"
tmp_dir=`echo $(mktemp -d)`

run_out='output.out'
EGH_out='EGH.out'

EGHrun_path=~/programs/EGHrun/

modes_fname='modes.dat'

time mpirun python $EGHrun_path/source/EGHrun.py -g $geom_xyz $geom_mol -rs $script_template -tdir $tmp_dir -run_out $run_out -EGH_out $EGH_out --calc_grad --calc_hess -modes $modes_fname -incr 0.01

exit 0
