#!/bin/bash
#SBATCH -n 32
#SBATCH -p gpus
#SBATCH -J alk2_pka_min
#SBATCH -o alk2_pka_min.o%j

export AMBERHOME=~/amber15mpi
source $AMBERHOME/amber.sh

module load openmpi/gcc/64/1.6.5

DO_MMPBSA="mpirun -n 32 $AMBERHOME/bin/pmemd.MPI -O "

struc_dir="structures"
output_dir="outfiles"
system_list=("R206H_vac" ) 

set -e
for ii in `seq 0 1 $((${#system_list[@]} - 1))`
do
	system=${system_list[$ii]}
	echo "minimizing system $system"
	$DO_MMPBSA -i min.mdin -p $struc_dir/$system.parm7 \
		-c $struc_dir/$system.rst7 -o $system.min.mdout \
		-r $outdir/$system.min.rst7 -ref $struc_dir/$system.rst7 \
		-cpin $system.cpin
done
