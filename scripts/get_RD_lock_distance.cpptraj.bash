#!/bin/bash
source /cm/shared/apps/amber18/amber.sh

system=$1
step=$2

topFile=${1}.parm7
trajFiles=${system}.${step}.[0-9][0-9].nc

maskStr=" :183@CG :204@CZ "

cpptraj $topFile << EOF
	trajin $trajFiles
	distance RD_Lock $maskStr out RD_Lock_Data.dat
	go
EOF
