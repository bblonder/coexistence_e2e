# Computation configs
# CORES = 4
# REPLICATES = 10
# GRID_POINTS = 20
# MIN_POINTS = 1e1
# MAX_POINTS = 1e4
# EXPERIMENTAL_DESIGNS = c("high-1", "high-2", "low-2", "low-3", "mixed", "prior")
CORES = 4
REPLICATES = 5
GRID_POINTS = 10
MIN_POINTS = 1e1
MAX_POINTS = 1e3
METHODS = c('naive','glv','rf','glv_rf_full')
# METHODS = c('rf', 'glv_rf_full')
EXPERIMENTAL_DESIGNS = c("mixed")
NUM_TEST = 100
MODEL_HYPERPARAMS = list(
  "num_factor_bins" = 99,
  "sequential_rf" = list(
    "iterations" = 3, 
    "nearest_k" = 10, 
    "bootstrap" = 10,
    "normalize" = TRUE,
    "score_weights" = list("uncertainty" = 1.0, "diversity" = 1.0, "density" = -1.0))
)

# Debugging override
if (DEBUG_MODE == TRUE) {
  CORES = 1
  REPLICATES = 2
  GRID_POINTS = 5
  MIN_POINTS = 1e1
  MAX_POINTS = 1e3
  EXPERIMENTAL_DESIGNS = c("mixed")
}