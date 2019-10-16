#!/bin/bash
#SBATCH -J alk_aa
#SBATCH -o alk_aa.o%j
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
struc_dir="."
output_dir="."

stepBase="prod"
step="prod_start"
for run_step in `seq -f "%02g" 0 1 4` 
do  prev_step=$step
	step=$stepBase.$run_step
	echo "--- --- ---"
	echo "Running $step for $system"
	pmemd.cuda -O -i $stepBase.mdin -c $output_dir/$system.$prev_step.rst7 -p $struc_dir/$system.parm7 \
		-o $system.$step.mdout  -r $output_dir/$system.$step.rst7 -x $output_dir/$system.$step.nc \
		-ref ${struc_dir}/${system}.${prev_step}.rst7 -inf ${output_dir}/${system}.${step}.mdinfo
done
echo "--- --- ---"
echo "Done!"
