#!/bin/bash

for i in {1..300}
do
sbatch --export arg1=500,arg2=60,arg3=6,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmdalt.cmd
sbatch --export arg1=500,arg2=60,arg3=2,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmdalt.cmd
sbatch --export arg1=500,arg2=60,arg3=3,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmdalt.cmd
sbatch --export arg1=500,arg2=60,arg3=4,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmdalt.cmd
sbatch --export arg1=500,arg2=60,arg3=5,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmdalt.cmd
done
