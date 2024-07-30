# Setup output directory
try(dir.create(file.path(getwd(), 'outputs/figures'), recursive = TRUE))
try(dir.create(file.path(getwd(), 'outputs/statistical'), recursive = TRUE))
directory_string = file.path(getwd(), 'outputs/statistical')

DEBUG_MODE = FALSE
source('src/configs.R')
source('src/coexistence_love.R')


# pick a dataset
cross_validate_environments <- function(assemblages_this, name, num_train_cases, num_train_env, num_bootstraps)
{
  num_species_this = assemblages_this %>% select(contains("action")) %>% ncol
  
  environment_this = assemblages_this %>% 
    select(contains("initial")) %>%
    apply(1, paste, collapse='*')
  
  num_train_max_this = min(table(environment_this))
  print(sprintf("num_train=%d num_train_max=%d num_env_max=%d",num_train_cases, num_train_max_this,length(table(environment_this))))
  
  # test in the environments that were not part of training
  results_this_dataset = rbindlist(lapply(1:num_bootstraps, function(i) {
    if (num_train_max_this < num_train_cases)
    {
      outcomes_novel_environments = data.frame(mae_test=NA)
      outcomes_training_environments = data.frame(mae_test=NA)
    }
    else
    {
      environment_train = sample(unique(environment_this), size=num_train_env)
      environment_test = setdiff(unique(environment_this), environment_train)
      
      print(sprintf('env_train=%s env_test=%s', paste(environment_train,collapse=","), paste(environment_test,collapse=",")))
      
      # eventually edit these to be by environment
      state_idxs_train_this = sample(which(environment_this %in% environment_train), size=num_train_cases)
      state_idxs_test_this = which(environment_this %in% environment_test)
      state_idxs_test_this_null = setdiff(which(environment_this %in% environment_train), state_idxs_train_this)
  
      if (length(state_idxs_test_this)>0)
      {
        outcomes_novel_environments = perform_prediction_experiment_single(predictor_variable="_abundance",
                                             assemblages=assemblages_this, 
                                             num_train=length(state_idxs_train_this),
                                             state_idxs_train=state_idxs_train_this,
                                             state_idxs_test=state_idxs_test_this,
                                             method='rf', 
                                             num_species=num_species_this)
      }
      else
      {
        outcomes_novel_environments = data.frame(mae_test=NA)
      }
      
      # test in the environments that were part of training
      outcomes_training_environments = perform_prediction_experiment_single(predictor_variable="_abundance",
                                           assemblages=assemblages_this, 
                                           num_train=length(state_idxs_train_this),
                                           state_idxs_train=state_idxs_train_this,
                                           state_idxs_test=state_idxs_test_this_null,
                                           method='rf', 
                                           num_species=num_species_this)
    }
    
    return(data.frame(name=name,
                      i=i,
                      num_train_cases=num_train_cases,
                      num_train_env=num_train_env,
                      mae_test_novel_environments=outcomes_novel_environments$mae_test,
                      mae_test_training_environments=outcomes_training_environments$mae_test))
  }))
  
  return(results_this_dataset)
}

cross_validate_environments_wrapper <- function(data, name, num_train_cases, num_train_env, num_bootstraps)
{
  params_this = expand.grid(num_train_cases=num_train_cases,  num_train_env=num_train_env, num_bootstraps = num_bootstraps)
  
  data_this = data %>% 
    mutate(state_idx=row_number()) # assume one row = one state
  
  results_this = rbindlist(pblapply(1:nrow(params_this), function(i) { 
    cross_validate_environments(assemblages_this = data_this,
                                name = name,
                                num_train_cases = params_this$num_train_cases[i],
                                num_train_env = params_this$num_train_env[i],
                                num_bootstraps = params_this$num_bootstraps[i])
    
    }))
  
  return(results_this)
}














## DO RUNS
NUM_TRAIN_CASES = generate_sample_size_sequences(MIN_POINTS, MAX_POINTS, GRID_POINTS, 200)

# ciliates
fn_ciliates = 'data/ciliates/data_ciliates.csv'
data_ciliates = read.csv(fn_ciliates)
results_ciliates = cross_validate_environments_wrapper(data = data_ciliates,
                                                       name='ciliates',
                                    num_train_cases = NUM_TRAIN_CASES,
                                    num_train_env = 1:6,
                                    num_bootstraps = 10)


# grasslands drought
fn_grassland_annual_plants_drought = 'data/grassland_annual_plants_drought/data_grassland_annual_plants_drought.csv'
data_grassland_annual_plants_drought = read.csv(fn_grassland_annual_plants_drought) %>%
  mutate(treatment.initial = factor(treatment.initial))
results_grassland_annual_plants_drought = cross_validate_environments_wrapper(data = data_grassland_annual_plants_drought,
                                                                              name='grassland_annual_plants_drought',
                                                                               num_train_cases = NUM_TRAIN_CASES,
                                                                               num_train_env = 1:2,
                                                                               num_bootstraps = 10)



fn_fruit_flies = 'data/fruit_flies/data_fruit_flies.csv'
data_fruit_flies = read.csv(fn_fruit_flies) %>%
  mutate(food.initial = factor(food.initial)) %>%
  mutate(temperature.initial = factor(temperature.initial))
data_fruit_flies[,1:28] = round(data_fruit_flies[,1:28]>0) # convert initial abundances to presence/absence

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

set.seed(1)
results_fruit_flies = cross_validate_environments_wrapper(data = data_fruit_flies,
                                                          name='fruit_flies',
                                                       num_train_cases = NUM_TRAIN_CASES,
                                                       num_train_env = 1:4,
                                                       num_bootstraps = 10)






source('src/dataset_stats.R')
df_stats = rbindlist(lapply(c(fn_ciliates, fn_grassland_annual_plants_drought, fn_fruit_flies), get_dataset_stats))


results_all = rbind(results_ciliates, 
                    results_grassland_annual_plants_drought, 
                    results_fruit_flies) %>% 
  pivot_longer(!name:num_train_env,names_to='test_mode') %>%
  left_join(df_stats,by='name') %>% 
  mutate(value_scaled = value / q_995)

write.csv(results_all, file='outputs/statistical/cross_validation.csv', row.names = FALSE)





# 
# results_ciliates  %>%
#   ggplot(aes(x=num_train_cases,y=value,color=test_mode)) + 
#   geom_boxplot() +
#   #scale_color_viridis_d() +
#   facet_grid(~num_train_env) +
#   theme_bw()
# 
# 
# results_grassland_annual_plants_drought %>% 
#   mutate(name=factor(name)) %>%
#   pivot_longer(!name:num_train_env,names_to='test_mode') %>%
#   ggplot(aes(x=num_train_cases,y=value,color=test_mode)) + 
#   geom_boxplot() +
#   #scale_color_viridis_d() +
#   facet_grid(~num_train_env) +
#   theme_bw()
# 
# results_fruit_flies %>% 
#   mutate(name=factor(name)) %>%
#   pivot_longer(!name:num_train_env,names_to='test_mode') %>%
#   ggplot(aes(x=num_train_cases,y=value,color=test_mode)) + 
#   geom_boxplot() +
#   #scale_color_viridis_d() +
#   facet_grid(~num_train_env) +
#   theme_bw()

