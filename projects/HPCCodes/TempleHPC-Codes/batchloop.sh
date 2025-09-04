#!/bin/sh
id_start=251
id_end=500
dir=/home/tuh36651/work/Spine6RyR/Runs
for ((x=id_start; x<=id_end; x++))
do
        qsub ${dir}/script_${x}.sh
done
