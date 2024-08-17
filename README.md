# coexistence_love
Community assembly outcome prediction and prioritization - Learning Outcomes Via Experiments (LOVE)

## Author
Benjamin Blonder (benjamin.blonder@berkeley.edu) and Michael H. Lim (michaelhlim@berkeley.edu)

## Datasets
All included datasets are obtained from public repositories and re-shared here. We do not claim ownership over any of these files. More information on data provenance is available in the Supporting Information of the accompanying manuscript.

## Getting started
Sequentially run the R scripts below in order to replicate all results from the main study.
* `1-predict.R`
* `2-plot.R`
* `3-plot_abundances.R`
* `4-cross_validate.R`
* `5-plot-cross-validate.R`

## Changing parameters
To change run parameters (e.g. `CORES`), see `src/config.R`. The default parameters for the first script will generate multiple gigabytes of outputs which are used to run the second script. Note that to run the first script with `CORES=20` requires >100 gigabytes of available memory due to the memory use required for some of the larger datasets.