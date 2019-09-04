#!/bin/bash
#SBATCH -J alk2_pka
#SBATCH -o alk2_pka.o%j
#SBATCH --partition=gpus
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --tasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --time=300:00:00

echo "starting up"

source /etc/profile.d/modules.sh
srun hostname -s | sort -u > slurm.hosts

module switch cuda55/toolkit/5.5.22 cuda70/toolkit/7.0.28

export CUDA_HOME=/cm/shared/apps/cuda70/toolkit/current

export AMBERHOME=/home/wesley/amber15c
source $AMBERHOME/amber.sh

cd /home/wesley/alk2_pKa/amber_structures
echo "entering directory: "`pwd`
export AMBER_JOBNAME=/home/wesley/alk2_pKa/amber_structures

system="R206H_vac"
prev_step="min"
struc_dir="structures"
output_dir="outfiles"

echo "running heating"
step="heat"
	pmemd.cuda -O -i $step.mdin -c $output_dir/$system.$prev_step.rst7 -p $struc_dir/$system.parm7 \
		-ref $struc_dir/$system.rst7 -cpin $system.cpin -o $system.$step.mdout \
		-r $output_dir/$system.$step.rst7 -x $output_dir/$system.$step.nc

echo "running equilibration"
prev_step=$step
step="equil"
	pmemd.cuda -O -i $step.mdin -c $output_dir/$system.$prev_step.rst7 -p $struc_dir/$system.parm7 \
		-cpin $system.cpin -o $system.$step.mdout -cpout $system.$step.cpout -cprestrt $system.$step.cpin \
		-r $output_dir/$system.$step.rst7 -x $output_dir/$system.$step.nc

