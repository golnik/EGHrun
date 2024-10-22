#!/bin/bash

#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 8
#SBATCH --time=00:30:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

gamess_bin="/home/ngolubev/Packages/software/gamess-uk/7.0/bin/gamess"

module add gcc/7.3.0

tdir="$TMPDIR/MP2_hessian/C3O2H2"		 # temporary root-directory

# normally no changes beyond this point

# define filenames
zmatrix=_zmatrix.dat
basisset=_basisset.dat

# go to temporary directory
mkdir -p $tdir
cp $zmatrix $tdir
cp $basisset $tdir
cd $tdir

# gamess-environment:
export ed3=dfile
export ed6=vfile

# some output
echo This calculation was done on `hostname`

#
#
#############################################################################
#
# SCF-Calculation:

echo "#begin<scf>"
$gamess_bin << EOF
title
$molec SCF calculation (basis $basis)
charge 0
`cat $zmatrix`
BASIS 6-31G
SCFTYPE MP2
RUNTYPE HESSIAN
EOF
echo "#end<scf>"

