library(ranger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(caret)
library(randomForestSRC)
library(parallel)
library(data.table)
library(e1071)
library(vegan)

# Set working directory
setwd("~/Dropbox/Berkeley/Research/Projects/Current/ecosys/coexistence_love")

# SETUP OUTPUT DIRECTORY
dir.create(file.path(getwd(), 'outputs/figures'), recursive = TRUE)
dir.create(file.path(getwd(), 'outputs/statistical'), recursive = TRUE)
directory_string = file.path(getwd(), 'outputs/statistical')

# Load helpers and settings
DEBUG_MODE = FALSE
source('src/configs.R')
source('src/coexistence_love.R')

# Perform analyses
set.seed(1)
data_soil_bacteria_8 = read.csv('data/friedman_gore/data_friedman_gore.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_soil_bacteria_8,
  dataset_name = 'soil_bacteria',
  num_species = 8, 
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 2) # this is an underestimate but should not cause problems

set.seed(1)
data_assemblages_M_11 = read.csv('data/glv/assemblages_M_11.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_assemblages_M_11,
  dataset_name = 'mouse_gut',
  num_species = 11,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)

set.seed(1)
data_annual_plant_18 = read.csv('data/annual_plant/assemblages_annual_plant_18.csv')
if (DEBUG_MODE==TRUE) {
 data_annual_plant_18 = data_annual_plant_18 %>% sample_n(2^14)
}
results = perform_prediction_experiment_full(
  directory_string,
  data_annual_plant_18,
  dataset_name = 'annual_plant',
  num_species = 18,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)

set.seed(1)
data_assemblages_H_12 = read.csv('data/glv/assemblages_H_12.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_assemblages_H_12,
  'human_gut',
  num_species = 12,
  method_list = METHODS,
  experimental_design_list = EXPERIMENTAL_DESIGNS,
  num_replicates_in_data = 1)

set.seed(1)
data_sortie_9_3 = read.csv('data/sortie/data_sortie.csv')
results = perform_prediction_experiment_full(
  directory_string,
  data_sortie_9_3,
  dataset_name = 'sortie-nd_plants',
  num_species = 9,
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
