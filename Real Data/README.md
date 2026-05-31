This folder contains the code required to reproduce the real data application presented in the paper *“Inferring Hub Nodes on Differential Score Matching Graphical Models”*.

The folder includes scripts for parameter tuning, model estimation, and visualization.

## Parameter Tuning

To reproduce the parameter tuning results, users should proceed as follows:

1. We find that the `lambda1` tuning results for normal conditional graphical model is unstable due to random data splitting. 
We obtain a robust tuning result by performing multiple data splits. 
To do this, submit the script `Tuning/tuning.R` by using `submit.sh` and `cmd.cmd` on a SLURM-based computing cluster to obtain the tuning results for normal conditional graphical model.

2. Run the `Tuning/NC_TunningSummary.R` script to obtain the final tuning results for `lambda1` and `lambda2` under normal conditional graphical model.

3. Run the `Tuning/Gaussian_and_Exp.R` script to obtain the tuning results for Gaussian and Exponential graphical models that are stable with respect to random data split.

## Model Estimation and Visualization

To directly estimate the model with given tuned parameters, users should proceed as follows:

1. Run the `RealData.R` script to obtain Figure 2 in the main manuscript and Figure F3 in the online supplementary material.