This folder contains the code required to reproduce all the simulation results presented in the paper *“Inferring Hub Nodes on Differential Score Matching Graphical Models”*.

Users are required to submit all the scripts in a SLURM-based computing cluster to collect the simulation results, and run the `Simu_Summary.R` script to obtain the corresponding tables and figures in the paper.

To submit the job, users are required to change line 6 and line 11 of the `cmd.cmd` script in each folder according to their own account and working directory.


## Repreducing Table 1

To reproduce Table 1 in the paper, users should proceed as follows:

1. Submit all the scripts in `NCHub_Null` and `ExpHub_Null` to the computing cluster.

2. Run the corresponding codes in `Simu_Summary.R`.

## Repreducing Figure 1

To reproduce Figure 1 in the paper, users should proceed as follows:

1. Submit all the scripts in `NCHub_Alt` and `ExpHub_Alt` to the computing cluster.

2. Run the corresponding codes in `Simu_Summary.R`.


## Repreducing Table 1 in the Online Supplementary Material

To reproduce Table 1 in the online supplementary material, users should proceed as follows:

1. Submit all the scripts in `NCSig` and `ExpSig` to the computing cluster.

2. Run the corresponding codes in `Simu_Summary.R`.




## Repreducing Tables 3--5 in the Online Supplementary Material

To reproduce Table 3--5 in the online supplementary material, users should proceed as follows:

1. Before running the Julia scripts, users should run the following commands in Julia:

```
(@v1.5) pkg> add https://github.com/mlakolar/KLIEPInference.jl
```

`using Pkg; Pkg.add.(["JLD", "Distributions", "ProximalBase", "CoordinateDescent", "StatsBase", "PyPlot", "IJulia", "Revise", "CSV", "DataFrames"])`

2. Submit all the scripts in `Gaussian_Compare` to the computing cluster.

3. Run the corresponding codes in `Simu_Summary.R`.

