# coexistence_e2e
End to end machine learning for coexistence outcome prediction

## Author
Benjamin Blonder (benjamin.blonder@berkeley.edu)

## Getting started
Sequentially run the R scripts below in order to replicate all results from the main study.
* `1 predict coexistence.R`* `2 plot coexistence.R`* `3 plot abundances.R`

Note that in the first file you may wish to change the run parameters (e.g. `CORES`) depending on the capabilities of your machine.

## Using the code on your own datasets
* Source the functions in `1 predict coexistence.R` and match the template function calls. `do_predictions()` is the main wrapper function. Dataset headers should be identical to the column headers used in the included datasets.