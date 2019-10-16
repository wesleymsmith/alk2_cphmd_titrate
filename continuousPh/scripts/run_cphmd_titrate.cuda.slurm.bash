#!/bin/bash
#SBATCH -J cphmd
#SBATCH -o cphmd.o%j
#SBATCH --partition=gpus
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --time=300:00:00

echo "starting up"

sysName=$1
rep=$2

repName=${sysName}.${rep}

numReps=20

source /etc/profile.d/modules.sh
srun hostname -s | sort -u > slurm.hosts

source /home/wesley/amber18/amber.sh

module load cuda91/toolkit/9.1.85
export CUDA_HOME=/cm/shared/apps/cuda91/toolkit/9.1.85

echo "running step 00 for $repName"
pmemd.cuda -O  -i $repName.in -c ${sysName}.rst7 -p ${sysName}.parm7 -phmdin cphmd.phmdin -phmdparm cphmd.phmdparm \
	-phmdout $repName.00.lambda -phmdrestrt $repName.00.phmdrst -o $repName.00.mdout -r $repName.00.rst -x $repName.00.nc \
	-inf $repName.00.mdinfo -l $repName.00.logfile
last_step=00
for ii in `seq -f "%02g" 1 1 4`
do	echo "running step $ii"
	pmemd.cuda -O  -i $repName.in -c $repName.${last_step}.rst -phmdstrt $repName.${last_step}.phmdrst \
	-p ${sysName}.parm7 -phmdin cphmd.phmdin -phmdparm cphmd.phmdparm \
	-phmdout $repName.${ii}.lambda -phmdrestrt $repName.${ii}.phmdrst -o $repName.${ii}.mdout -r $repName.${ii}.rst -x $repName.${ii}.nc \
	-inf $repName.${ii}.mdinfo -l $repName.${ii}.logfile
	last_step=$ii
done

echo "done"
