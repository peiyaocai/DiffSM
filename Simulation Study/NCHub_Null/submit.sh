#!/bin/bash

for i in {1..300}
do



sbatch --export arg1=400,arg2=50,arg3=1,arg4=$i,arg5=5 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=500,arg2=60,arg3=1,arg4=$i,arg5=5 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=600,arg2=70,arg3=1,arg4=$i,arg5=5 -o /dev/null -e /dev/null cmd.cmd




done
