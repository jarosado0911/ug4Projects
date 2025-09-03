#!/bin/sh 
#PBS -l walltime=00:60:00 
#PBS -N Spine5_activeRyR
#PBS -q normal 
#PBS -l nodes=2:ppn=28 
#PBS -e error.txt 
#PBS -o output.txt 
cd $PBS_O_WORKDIR 

mpirun -np 56 ugshell -ex /home/tuh36651/work/Spine5_activeER/reconstructed_spine_wER.lua -grid /home/tuh36651/work/Spine5_activeER/Spine5_wER.ugx -numRefs 0 -caInflux 0.0086 -tstep 5e-6 -endTime 0.015 -outName /home/tuh36651/work/Spine5_activeER/output -solver GS -setting ryr -ryrDensity 4.0 -vtk -pstep 0.0001 -minDef 1e-12 -freq 1
