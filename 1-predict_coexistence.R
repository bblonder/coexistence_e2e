# Setup output directory
setwd("~/Documents/coexistence_love")
dir.create(file.path(getwd(), 'outputs/figures'), recursive = TRUE)
dir.create(file.path(getwd(), 'outputs/statistical_new'), recursive = TRUE)
directory_string = file.path(getwd(), 'outputs/statistical_new')

# Load helpers and settings
DEBUG_MODE = TRUE
source('src/configs.R')
source('src/coexistence_love.R')

# if on cluster
# CORES <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))

# Perform analyses
set.seed(1)
data_soil_bacteria_8 = read.csv('data/soil_bacteria/data_soil_bacteria.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_soil_bacteria_8,
  dataset_name = 'soil_bacteria',
  num_species = 8, 
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 2) # this is an underestimate but should not cause problems

set.seed(1)
data_mouse_gut_11 = read.csv('data/human_and_mouse_gut/data_mouse_gut.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_mouse_gut_11,
  dataset_name = 'mouse_gut',
  num_species = 11,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)

set.seed(1)
data_human_gut_12 = read.csv('data/human_and_mouse_gut/data_human_gut.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_human_gut_12,
  'human_gut',
  num_species = 12,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)

set.seed(1)
data_ciliates = read.csv('data/ciliates/data_ciliates.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_ciliates,
  dataset_name = 'ciliates',
  num_species = 6,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 3)

set.seed(1)
data_fly_5 = read.csv('data/fly/data_fly.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_fly_5,
  dataset_name = 'fly_gut',
  num_species = 5,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 48)

set.seed(1)
data_assemblages_cedar_creek_18 = read.csv('data/cedar_creek/cedar_creek_2018.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_assemblages_cedar_creek_18,
  dataset_name = 'cedar_creek_plants',
  num_species = 18,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)

set.seed(1)
data_annual_plant_18 = read.csv('data/annual_plant/assemblages_annual_plant_18.csv')
if (DEBUG_MODE==TRUE) {
  data_annual_plant_18 = data_annual_plant_18 %>% sample_n(2^14)
}
METHODS = setdiff(METHODS, 'sequential_rf') # no sequential RF due to large parameter space
results = perform_prediction_experiment_full(
  directory_string,
  data_annual_plant_18,
  dataset_name = 'annual_plant',
  num_species = 18,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)
