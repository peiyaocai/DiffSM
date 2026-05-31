#!/bin/bash

for i in {1..15}
do
sbatch --export arg1=400,arg2=50,arg3=1,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=400,arg2=50,arg3=2,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=400,arg2=50,arg3=3,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=400,arg2=50,arg3=4,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=400,arg2=50,arg3=5,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd

sbatch --export arg1=500,arg2=50,arg3=1,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=500,arg2=50,arg3=2,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=500,arg2=50,arg3=3,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=500,arg2=50,arg3=4,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=500,arg2=50,arg3=5,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd

sbatch --export arg1=600,arg2=50,arg3=1,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=600,arg2=50,arg3=2,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=600,arg2=50,arg3=3,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=600,arg2=50,arg3=4,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd
sbatch --export arg1=600,arg2=50,arg3=5,arg4=$i,arg5=6 -o /dev/null -e /dev/null cmd.cmd

done
