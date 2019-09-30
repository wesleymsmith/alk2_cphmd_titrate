#!/bin/bash
#SBATCH -J cphmdremd
#SBATCH -o cphmdremd.o%j
#SBATCH --partition=gpus
#SBATCH --get-user-env
#SBATCH --nodes=5
#SBATCH --tasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --time=300:00:00

sysName=$1
runNum=$2

echo "starting up"

numReps=`cat groupfile.template | wc -l`

source /etc/profile.d/modules.sh
srun hostname -s | sort -u > slurm.hosts

source /home/wesley/amber18/amber.sh

module load cuda91/toolkit/9.1.85

if [ "$runNum" -eq "00" ]
then
	cat groupfile.template | \
		sed "s/rep[0-9][0-9].__RUN0__.rst/${sysName}.min.rst7/g" | \
		sed "s/-phmdstrt rep[0-9][0-9].__RUN0__.phmdrst//g" | \
		sed "s/__SYSTEM__/$sysName/g" | \
		sed "s/__RUN1__/00/g" > groupfile.$runNum
else
	lastRun=`printf "%02g" $(($runNum - 1))`
	cat groupfile.template | \
		sed "s/__RUN0__/$lastRun/g" | \
		sed "s/__RUN1__/$runNum/g" | \
		sed "s/__SYSTEM__/$sysName/g" > groupfile.$runNum
fi

if [ "$3" == "dry_run" ]
then
	echo "here is the groupfile:"
	echo "--- --- ---"
	cat groupfile.$runNum
	echo "--- --- ---"
else
	mpirun -np $numReps pmemd.cuda.MPI -ng $numReps -groupfile groupfile.$runNum -rem 4
fi
echo "done"
