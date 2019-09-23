# Biased Pre-Processing Strategies for Imbalanced Spatio-Temporal Forecasting

This repository contains the research compendium of the conference article:

Oliveira M, Moniz N, Torgo L, Santos Costa V (2019). “Biased Resampling Strategies for Imbalanced Spatio-Temporal Forecasting.” In _Proceedings of the IEEE International Conference on Data Science and Advanced Analytics (DSAA)_.

You are free to use and/or adapt the code we freely provide. However, we do require that if you do that you cite the paper where these results and code were published.

If you adapt the code to your own needs, you are also required to maintain information on your code concerning the original source of the code (e.g. the URL of this repository) and a reference to the original paper.

## Prerequisites

To install this package, run:

```
library(devtools)  # You need to install this package!
install_github("mrfoliveira/STResampling-DSAA2019",ref="master")
```

If this previous install_github calls somehow fail (there are reports of problems with different libcurl library version on Linux hosts) you may try in alternative the following in R:

```
library(devtools)
install_git("mrfoliveira/STResampling-DSAA2019",ref="master")
```

## Reproducing experiments

To obtain all results shown in the article, run the following lines from the package directory:

```
library(STResamplingDSAA)
source("inst/generate_inds.R")
source("inst/exps_internalTuning.R")
source("inst/exps_externalPrequential.R")
```

## Contents

The repository is organized as follows:

* **inst/** - Contains scripts for spatio-temporal data, scripts for launching experiments, summarized experimental results, and a report generating all the figures and tables in the conference article
* man/ - Contains function documentation
* **R/** - Contains R package code
* STResamplingDSAA_0.5-1.pdf - reference manual

#### inst
* CITATION - file containing citation information
* exps_externalPrequential.R - script to launch experiments using a prequential evaluation method
* exps_internalTuning.R - script to launch experiments with internal parameter tuning
* generate_inds.R - script to transform data sets and generate spatio-temporal indicators
* notebook.Rmd - Rmarkdown file capable of generating a HTML report containing figures and tables present in the paper
* sumRes_externalPrequential.Rdata - object containing results of experiments with all sets of parameters
* sumRes_internalTuning.Rdata - object containing results of experiments with internally tuned parameters
* extdata/dfs.Rdata - object containing list of data sets used
* figs/ - directory containing all figures shown in the paper and presentation

#### R

* eval_framework.R - functions needed for experimental evalution: error estimation methods and
error metric calculation functions, depends on utils.R
* sampling_methods - resampling functions with a spatio-temporal bias
* sampling_weight.R - functions to calculate biased re-sampling weights
* st_indicators.R - functions to calculate spatio-temporal indicators
* UBL_sampling_methods.R - resampling functions almost identical to the functions in R package UBL (slight difference in the way one argument functions)
* utils.R - general utility functions
* workflows.R - functions that define a learning/predicting workflow
