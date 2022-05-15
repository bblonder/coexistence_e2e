# coexistence_love
Machine learning for coexistence outcome prediction and prioritization

## Author
Benjamin Blonder (benjamin.blonder@berkeley.edu)

## Datasets
All included datasets are obtained from public repositories and re-shared here. I do not claim ownership over any of these files. More information on data provenance is available in the Supporting Information of the accompanying manuscript.

## Getting started
Sequentially run the R scripts below in order to replicate all results from the main study.
* `1 predict coexistence.R`
* `2 plot coexistence.R`
* `3 plot abundances.R`

Note that in the first file you may wish to change the run parameters (e.g. `CORES`) depending on the capabilities of your machine. The default parameters for the first script will generate approximately 6GB of outputs which are used to run the second script.