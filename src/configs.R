# Computation configs
CORES = 4
REPLICATES = 10
GRID_POINTS = 20
MIN_POINTS = 1e1
MAX_POINTS = 1e4
METHODS = c('naive','rf')
EXPERIMENTAL_DESIGNS = c("high-1","high-2", "low-2","low-3", "mixed")
NUM_TEST = 100

# Debugging override
if (DEBUG_MODE == TRUE) {
  CORES = 1
  REPLICATES = 1
  GRID_POINTS = 1
  MIN_POINTS = 500
  MAX_POINTS = 500
  METHODS = c('rf')
  EXPERIMENTAL_DESIGNS = c("mixed")
}