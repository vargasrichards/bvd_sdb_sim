## Simulation of Safe and Dignified Burial for control of Bundibugyo Ebola virus epidemics
### Description of input datasets and R analysis scripts
May 2026

## General description
This repository contains data and R scripts needed to replicate the above analysis. The input datasets required are found in the `\in` folder, and are read automatically when the code is run.
All the code is found in the `\code` folder. To replicate the analysis, follow these steps:
* Download and unzip the repository to any folder in your computer (other than the Downloads folder, which usually gets wiped automatically). The folder is identified automatically when the code is run.
* Download R and RStudio (see download links on [https://posit.co/download/rstudio-desktop/]). While R is sufficient to run the analysis, it is recommended to instead run the scripts from the RStudio interface.
* Open and run the entire `00_control_script.R` script (just press Alt+Ctrl+R). This will create an `\out` folder with further sub-folders, to which output tables and graphs will be saved automatically. As this scripts calls all the others, it alone is sufficient to replicate the analysis, but note below steps if you wish to alter the analysis.

## Description of input files
* `evd_sdb_sim_parameters.xlsx` contains a list of model parameters which that user can modify as wished. Note that the `n_sim` parameter heavily affects computing power requirements: on a cheap 2020 Laptop, 1000 simulation runs are OK for about 20-30 individual scenarios (and take <10 min to run), but anything beyond this may require a powerful PC or a high-performance cluster.
* `out_dose_resp_rn_p_success_hi.csv` and `out_dose_resp_rn_p_success_pw.csv` are generated from a separate analysis, described in Checchi et al. (2025) (https://researchonline.lshtm.ac.uk/id/eprint/4678756/) and which the user can reproduce entirely from the following repository: https://github.com/francescochecchi/evd_drc_sdb_effect. The two files will appear in the `/out` directory created after running that analysis. Each file contains the model-predicted net reproduction number of EVD at different levels of Safe and Dignified Burial (SDB) effectiveness (proportion of successful instances or <p_success> column), according to the Hirano-Imbens method (file name ending with `_hi.csv`) or the propensity weights (`_pw.csv`) method. The other columns report the mean reduction in R and the lower and 95% confidence interval estimates.

## Description of R scripts
* `00_control_script.R` loads necessary packages, sets a few general parameters and calls other scripts.
* `01_prepare_model.R`loads model parameters, sets scenario values and initialises a stochastic transmission dynamic model using the excellent <pomp> package (https://kingaa.github.io/pomp/).
* `02_main_analysis.R` carries out the main analysis. This includes running the main scenarios to be highlighted and visualising/saving some results.
* `03_sens_analyses.R` carries out two sensitivity analyses: (i) effect of varying the number of starting (seed) cases; and (ii) two extreme scenarios (SDB vs. no SDB). 
