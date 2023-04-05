# Computation configs
# CORES = 4
# REPLICATES = 10
# GRID_POINTS = 20
# MIN_POINTS = 1e1
# MAX_POINTS = 1e4
CORES = 4
REPLICATES = 10
GRID_POINTS = 20
MIN_POINTS = 1e1
MAX_POINTS = 1e4
# METHODS = c('rf')
# EXPERIMENTAL_DESIGNS = c("high-1","high-2", "low-2","low-3", "mixed", "prior")
METHODS = c('naive','rf','sequential_rf')
EXPERIMENTAL_DESIGNS = c("mixed", "prior")
NUM_TEST = 100
MODEL_HYPERPARAMS = list(
  "sequential_rf" = list(
    "iterations" = 3, 
    "nearest_k" = 10, 
    "bootstrap" = 10,
    "normalize" = TRUE,
    "score_weights" = list("uncertainty" = 1.0, "diversity" = 1.0, "density" = -1.0))
)

# Debugging override
if (DEBUG_MODE == TRUE) {
  CORES = 4
  REPLICATES = 1
  GRID_POINTS = 5
  MIN_POINTS = 1e1
  MAX_POINTS = 1e4
  # METHODS = c('naive','rf','sequential_rf')
  METHODS = c('sequential_rf')
  EXPERIMENTAL_DESIGNS = c("prior")
}