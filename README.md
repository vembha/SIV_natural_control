# **Code repository for the publication:**
# Antiviral capacity of the early CD8 T-cell response is predictive of natural control of SIV infection: Learning _in vivo_ dynamics using _ex vivo_ data
### URL: [bioRxiv version](https://www.biorxiv.org/content/10.1101/2023.10.13.562306v1)

#### Numerical packages and languages used: Julia 1.7.3; Monolix 2021R1
<br/>
<br/>

This repository contains three folders named `figure_codes`, `figure_data_files`, and `monolix_files` which can be used to generate the figures in the manuscript. The `julia_packages.txt` file lists all the packages used in Julia 1.7.3 to perform the calculations. A description of the files in the individual folders is given below:

#### 1. `monolix_files`
These files fit the best-fit model of our study to the longitudinal data and infer the parameter values. The file `masterfile_model_1.mlxtran` is the Monolix file that should be run to obtain these estimates. While `modelFile_model_1.txt` codes the model, the `master_data_RNA_DNA_p27_kE0.csv` is the raw longitudinal data file; both files are formatted as per the requirement of Monolix. Finally, the `nohup.txt` is the output obtained by running the `.mlxtran` file.

#### 2. `figure_data_files`
This folder contains all the `.xlsx` and `.csv` files required to plot the figures. The names of the files are self-explanatory: files with the suffix `parameters_` have the parameter values estimated for different models employed in the study; the `master_data_` suffix refers to files with the raw longitudinal data--both _in vivo_ and _ex vivo_ models; and the `modified_DNA_post_peak.xlsx` is the SIV DNA data without the measurements before the peak in the data.

#### 3. `figure_codes`
This folder contains the Julia files that can be run to reproduce the figures; names are self-explanatory.

#### 4. `julia_packages
This text file provides the versions of the Julia packages used for this project.
