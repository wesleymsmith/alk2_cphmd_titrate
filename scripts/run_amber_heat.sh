#!/bin/bash
#SBATCH -J cph_heat
#SBATCH -o cph_heat.o%j
#SBATCH --partition=gpus
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=300:00:00

echo "starting up"

source /etc/profile.d/modules.sh
srun hostname -s | sort -u > slurm.hosts

module switch cuda55/toolkit/5.5.22 cuda91/toolkit/9.1.85

#export CUDA_HOME=/cm/shared/apps/cuda70/toolkit/current

export AMBERHOME=/home/wesley/amber18
source $AMBERHOME/amber.sh

#cd /home/wesley/alk2_pKa/amber_structures
#echo "entering directory: "`pwd`
#export AMBER_JOBNAME=/home/wesley/alk2_pKa/amber_structures

system=$1
step="min"
struc_dir=$2
output_dir=$3

echo "running heating"
prev_step=$step
step="heat"
	pmemd.cuda -O -i $step.mdin -c $output_dir/$system.$prev_step.rst7 -p $struc_dir/$system.parm7 \
		-cpin $system.cpin -o $system.$step.mdout -cpout $system.$step.cpout -cprestrt $system.$step.cpin \
		-r $output_dir/$system.$step.rst7 -x $output_dir/$system.$step.nc \
		-ref ${struc_dir}/${system}.${prev_step}.rst7 -inf ${output_dir}/${system}.${step}.mdinfo

