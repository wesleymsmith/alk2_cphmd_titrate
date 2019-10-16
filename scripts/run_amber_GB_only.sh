#!/bin/bash
#SBATCH -J alk_gb
#SBATCH -o alk_gb.o%j
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

export AMBERHOME=/cm/shared/apps/amber18
source $AMBERHOME/amber.sh

#cd /home/wesley/alk2_pKa/amber_structures
#echo "entering directory: "`pwd`
#export AMBER_JOBNAME=/home/wesley/alk2_pKa/amber_structures

system=$1
window_num=$2

output_dir="."
struc_dir="."

echo "running equilibration"

if [ "$window_num" -eq "00" ]
then
	prev_step="start"
else
	prev_step=`printf "%02g" $(($window_num -1))`
fi
step="prod.${window_num}"
	pmemd.cuda -O -i $system.mdin -c $output_dir/$system.$prev_step.rst7 -p $struc_dir/$system.parm7 \
		 -o $system.$step.mdout -r $output_dir/$system.$step.rst7 -x $output_dir/$system.$step.nc \
		-ref ${struc_dir}/${system}.${prev_step}.rst7 -inf ${output_dir}/${system}.${step}.mdinfo

