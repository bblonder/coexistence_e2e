# Setup output directory
try(dir.create(file.path(getwd(), 'outputs/figures'), recursive = TRUE))
try(dir.create(file.path(getwd(), 'outputs/statistical'), recursive = TRUE))
directory_string = file.path(getwd(), 'outputs/statistical')

# Load helpers and settings
DEBUG_MODE = FALSE
source('src/configs.R')
source('src/coexistence_love.R')

# if on cluster
# CORES <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))

# Perform analyses



### MULTIPLE ENVIRONMENTS
set.seed(1)
data_fruit_flies = read.csv('data/fruit_flies/data_fruit_flies.csv') %>%
  mutate(food.initial = factor(food.initial)) %>%
  mutate(temperature.initial = factor(temperature.initial))

data_fruit_flies[,1:28] = round(data_fruit_flies[,1:28]>0) # convert initial abundances to presence/absence

# drop 7 lowest-abundance in outcome species
species_lowest_abundance = data_fruit_flies %>% 
  select(contains("outcome")) %>% 
  colMeans %>% 
  sort %>% 
  head(7) %>%
  names %>%
  gsub("\\.outcome","",.)

data_fruit_flies = data_fruit_flies %>%
  select(!contains(species_lowest_abundance)) %>%
  as.data.frame

results = perform_prediction_experiment_full(
  directory_string,
  data_fruit_flies,
  dataset_name = 'fruit_flies',
  num_species = length(grep("outcome",names(data_fruit_flies))), # in full dataset, should be 28
  method_list = c('rf','sequential_rf','naive'),
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1) # in full dataset, this is maximum # of replicates, reflects the last 30 rows of the data file where there are varying abundances (which are ignored by this code run)




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
data_grassland_annual_plants_drought = read.csv('data/grassland_annual_plants_drought/data_grassland_annual_plants_drought.csv') %>%
  mutate(treatment.initial = factor(treatment.initial))
results = perform_prediction_experiment_full(
  directory_string,
  data_grassland_annual_plants_drought,
  dataset_name = 'grassland_annual_plants_drought',
  num_species = 6,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)





# SINGLE ENVIRONMENTS

set.seed(1)
data_soil_bacteria = read.csv('data/soil_bacteria/data_soil_bacteria.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_soil_bacteria,
  dataset_name = 'soil_bacteria',
  num_species = 8, 
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 2) # this is an underestimate but should not cause problems

set.seed(1)
data_mouse_gut = read.csv('data/human_and_mouse_gut/data_mouse_gut.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_mouse_gut,
  dataset_name = 'mouse_gut',
  num_species = 11,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)

set.seed(1)
data_human_gut = read.csv('data/human_and_mouse_gut/data_human_gut.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_human_gut,
  'human_gut',
  num_species = 12,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)


set.seed(1)
data_fly_gut = read.csv('data/fly_gut/data_fly_gut.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_fly_gut,
  dataset_name = 'fly_gut',
  num_species = 5,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 48)

set.seed(1)
data_prairie_plants = read.csv('data/prairie_plants/data_prairie_plants.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_prairie_plants,
  dataset_name = 'prairie_plants',
  num_species = 18,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)

set.seed(1)
data_forest_trees = read.csv('data/forest_trees/data_forest_trees.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_forest_trees,
  dataset_name = 'forest_trees',
  num_species = 9,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 3)


# this does not run on my cluster - too much memory to run in parallel. had to put on HPC. 
# run carefully on a local machine or reduce # of cores in src/configs.R ...
set.seed(1)
data_grassland_annual_plants = read.csv('data/grassland_annual_plants/data_grassland_annual_plants.csv')
if (DEBUG_MODE==TRUE) {
  data_grassland_annual_plants = data_grassland_annual_plants %>% sample_n(2^14)
}
results = perform_prediction_experiment_full(
  directory_string,
  data_grassland_annual_plants,
  dataset_name = 'grassland_annual_plants',
  num_species = 18,
  method_list = setdiff(METHODS, 'sequential_rf'), # no sequential RF due to large parameter space
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)


