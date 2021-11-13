# Paper Repository

This repository contains the code to reproduce the numerical results in our paper - Sensitivity Analysis of Individual Treatment Effects: A Robust Conformal Inference Approach.

The following R packages are required to be installed: [grf](https://grf-labs.github.io/grf/), [tidyverse](https://www.tidyverse.org/).





## Folders

- `simulations/`: R codes for simulations in Section 6. 

- `realdata/`: R codes for real data analysis in Section 7. 
- `bash/`: bash files run the simualtions in batch mode. 
- `utils/`: basic functions for simulations and real data analysis.



## Running simulations

Each file in `simulations/` folder implements one run of the simulation. 

### Single run

##### Counterfactual prediction (Sec 6.1)

Files starting with `pred_` are for counterfactual prediction in Section 6.1. Each of them take the inputs `--p` for the dimension, `--n` for the sample size, `--alpha_id` for the coverage target, `--gamma_id` for the confounding level and `--seed` for the random seed. 

Simulations in the paper take dimension 4 and 20, sample size 500, 2000 and 5000, `alpha_id` from 1 to 9 corresponding to coverage from 0.9 to 0.1, `gamma_id` from 1 to 5 corresponding to confounding level 1.5, 2, 2.5, 3, 5, and seed from 1 to 1000. The result of each run will be stored in (created) `results/` folder with file name indicating the algorithm employed and all these configurations.

For example, to execute a single run of counterfactual prediction using the marginally valid procedure with ground truth of bounds, covariate dimension 4, sample size 2000, coverage 0.9, confounding level 1.5 and random seed 5, one can run the following script:

```
cd simulations
Rscript pred_mgn.R 4 2000 1 1 5 
```

The marginally valid procedure with estimated bounds and same configurations can be executed by:

```
Rscript pred_mgn_est.R 4 2000 1 1 5 
```

The PAC-type procedure with ground truth of bounds / estimated bounds and same configurations can be executed by:

```
Rscript pred_pac.R 4 2000 1 1 5 
Rscript pred_pac_est.R 4 2000 1 1 5 
```

It stores the summary of experiment (coverage, estimation error and configurations) in `results/simulation/` folder, with file name showing all information of configurations.

##### Sensitivity analysis (Sec 6.2)

Files starting with `sens_` are for sensitivity analysis in Section 6.2.  Each of them take the inputs `--gamma_id` for the confounding level, `--diff_id` for the effect size, `ite` indicating random or fixed ITE, and `--seed` for the random seed. 

Simulations in the paper take `gamma_id` from 1 to 6 corresponding to confounding level 1, 1.2, 1.4 ,..., 2, `diff_id` from 1 to 5 corresponding to effect size -1, -0.5, 0, 0.5, 1, and seed from 1to 1000. The result of each run will be stored in `results/` folder with file name indicating the algorithm employed and all these configurations.

For example, to execute a single run of sensitivity analysis using the marginally valid procedure with ground truth of bounds, confounding level 1.4, effect size 0.5, random ITE and random seed 5, one can run the following script:

```
cd simulations
Rscript sens_mgn.R 3 4 1 5
```

It stores the summary of this single run (summary of coverage starting with `sens_mgn_res_`, ground truth of ITE starting with `sens_mgn_ITE_`, obtained Gamma-values on test samples starting with `sens_mgn_gvalue_`) in `results/simulation/` folder with all information of configurations.

One run using the marginally valid procedure with estimated bounds and the same configuration can be executed by:

```
Rscript sens_mgn_est.R 3 4 1 5
```

The PAC-type procedure with ground truth or estimated bounds and the same configuration except for fixed ITE can be executed by:

```
Rscript sens_pac.R 3 4 2 5
Rscript sens_pac_est.R 3 4 2 5
```

It stores the summary of this single run similarly in `results/simulation/` folder with all information of configurations.

### Batch submission 

The simulations can also be submitted in a batch mode on computing clusters, using bash files in `bash/` folder (might need modification according to the configurations of computing clusters). We implement all simulation configurations with seed from 1 to 100 (experiments in the paper use 1 to 1000). 

For example, to submit the jobs for counterfactual prediction with marginally valid procedure (both with ground truth and estimated bounds), direct to `bash/` folder and run 

```
sh run_simu_pred_mgn.sh
```

The random seeds can be changed in the bash files. It stores the outputs in separate files, in the same way for each single run. 

## Running real data analysis

Each file in `realdata/` folder implements one run of the real data analysis. 

### Single run

##### Counterfactual prediction (Sec 7.1)

Files starting with `syn_` are for counterfactual prediction of semi-real data in Section 7.1. We implement the two procedures with estimated bounds. 

To run the experiments, first direct to `realdata/` folder and generate the synthetic population with

```
cd realdata
Rscript acic_pred_generate.R
```

It will store an `.RData` file for subsequent experiments. 

The procedures take inputs `alpha_id` from 1 to 9 corresponding to target coverage from 0.1 to 0.9, `gamma_id` from 1 to 5 corresponding to true confounding level 1.5, 2, 2.5, 3, 5, and `seed` for random seed.

For example, to execute a single run of counterfactual prediction on semi-real data with marginally valid procedure, coverage 0.8, confounding level 2 and random seed 3, run the following command:

```
Rscript syn_pred_mgn.R 2 2 3
```

It stores the summary of results (coverage, average length of predictive interval, etc) in `results/realdata/` folder, with file name `syn_pred_marginal_` followed by all relevant configurations. 

The following command runs the PAC procedure with the same configuration:

```
Rscript syn_pred_pac.R 2 2 3
```

It stores the summary of results (coverage, average length of predictive interval, etc) in `results/realdata/` folder, with file name `syn_pred_pac_` followed by all relevant configurations. 

##### Sensitivity analysis

Files starting with `sens_` are for sensitivity analysis on real data in Section 7.2. Each of them takes inputs `pos` indicating whether we test for positive ITE (`pos=1`, null hypothesis: ITE <= 0 and Gamma* <= Gamma) or negative ITE (`pos=2`), `alpha_id` for target coverage and `seed` for random seed (for sample splitting). 

For example, to test for positive ITE with marginally-valid procedure at confidence level 0.8 and random seed 10, run

```
cd realdata
Rscript sens_mgn.R 1 2 10
```

It stores the covariate, treatment and responses of test data along with the obtained Gamma-values in `/results/realdata/` folder. The file name is `sens_positive_mgn_` followed by information of configurations. 

The same configuration with PAC-type procedure can be executed by 

```
Rscript sens_pac.R 1 2 10
```

It stores the covariate, treatment and responses of test data along with the obtained Gamma-values in `/results/realdata/` folder. The file name is `sens_positive_pac_` followed by information of configurations. 

### Batch submission 

Experiments on real data can be submitted in batch mode using bash files in `bash/` folder (might need modification according to the configurations of computing clusters). We implement all simulation configurations with seed from 1 to 100, which can be modified as well (the paper uses 1 to 1000). 

For example, to run all configurations using marginally valid procedure with random seed from 1 to 100, one can direct to `bash/` folder and run the following command:

```
sh run_real_pred_mgn.sh
```

To run all configurations using PAC-type procedure with random seed from 1 to 100, one can direct to `bash/` folder and run the following command:

```
sh run_real_pred_pac.sh
```

It stores outputs of all single runs in `results/realdata/` in the same way. 

For sensitivity analysis using marginally valid procedure with all configurations of confidence level and hypothesis type (10 random splits), direct to `bash/` folder and run the command:

```
sh run_real_sens_mgn.sh
```

For sensitivity analysis using PAC-type procedures, direct to `bash/` folder and run the command:

```
sh run_real_sens_pac.sh
```

Outputs of all single runs are stored in `results/realdata/` in the same way. 