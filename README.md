# Biased Pre-Processing Strategies for Imbalanced Spatio-Temporal Forecasting

This repository contains the research compendium of the conference article:

**Oliveira, M., Moniz, N., Torgo, L., & Costa, V. S. (2019, October). Biased Resampling Strategies for Imbalanced Spatio-Temporal Forecasting. In _2019 IEEE International Conference on Data Science and Advanced Analytics (DSAA)_ (pp. 100-109). IEEE. doi: [10.1109/DSAA.2019.00024](https://doi.org/10.1109/DSAA.2019.00024)**

BibTeX citation:

@inproceedings{Oliveira2019, author = {Mariana Oliveira and Nuno Moniz and Lu{\'{\i}}s Torgo and V{\'{\i}}tor Santos Costa}, editor = {Lisa Singh and Richard D. De Veaux and George Karypis and Francesco Bonchi and Jennifer Hill}, title = {Biased Resampling Strategies for Imbalanced Spatio-Temporal Forecasting}, booktitle = {2019 {IEEE} International Conference on Data Science and Advanced Analytics, {DSAA} 2019, Washington, DC, USA, October 5-8, 2019}, pages = {100--109}, publisher = {{IEEE}}, year = {2019}, url = {https://doi.org/10.1109/DSAA.2019.00024}, doi = {10.1109/DSAA.2019.00024}}

You are free to use and/or adapt the code we freely provide. However, we do require that if you do that you cite the paper where these results and code were published.

If you adapt the code to your own needs, you are also required to maintain information on your code concerning the original source of the code (e.g. the URL of this repository) and a reference to the original paper.

Other supplementary material (e.g., slides of conference presentation) available at https://www.dcc.fc.up.pt/~moliveira/publication/19-dsaa-biased-resampling-spatiotemporal/.

##### **Notice**

In Section IV.B2 of the article, we report that the sets of parameters tested in Section V are _u_ &isin; \{.2, .4, .6, .8, .95\} and _o_ &isin; \{.5, 1, 2, 3, 4\}. In fact, the results reported in Sections V.A and V.B for _Internally Tuned Parameters_ were obtained with _u_ &isin; \{.1, .2, .3, .4, .5, .6, .7, .8, .9\} and _o_ &isin; \{.1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 2\}. The results of running the internal tuning experiments with the initially reported sets of parameters does not alter the conclusions of our work.

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
PATH <- system.file("inst/", package="STResamplingDSAA")
source(paste0(PATH, "/generate_inds.R"))
source(paste0(PATH, "/exps_internalTuning.R"))
source(paste0(PATH, "/exps_externalPrequential.R"))
```

To generate an HTML report containing all figures and tables in the article, run:

```
library(STResamplingDSAA)
knitr::knit(system.file("inst/notebook.Rmd", package="STResamplingDSAA"))
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
