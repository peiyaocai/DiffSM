#!/bin/bash

for i in {1..300}
do
sbatch --export arg1=$i -o /dev/null -e /dev/null cmd.cmd
done
