#!/bin/bash -l

#SBATCH --nodes 1
#SBATCH --ntasks-per-node 10
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

conda info
which python

which mpirun

cat << EOF > test.py
from mpi4py import MPI

# Get the default communicator, which includes all processes
comm = MPI.COMM_WORLD

# Get the rank of the current process (its unique ID)
rank = comm.Get_rank()

# Get the total number of processes in the communicator
size = comm.Get_size()

# Get the name of the processor where the current process is running
name = MPI.Get_processor_name()

# Print a message from each process
print(f"Hello World! I am process {rank} of {size} on {name}.")
EOF

mpirun python test.py

rm -f test.py

exit 0
