#!/bin/sh
id_start=0
id_end=1
dir=/home/tuh36651/ug4/runs/Spine11ryrNeck1/Runs
for ((x=id_start; x<=id_end; x++))
do
        sbatch ${dir}/script_${x}.sb
done
