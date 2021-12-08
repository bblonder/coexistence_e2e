library(ranger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(caret)
library(randomForestSRC)
library(parallel)
library(data.table)
library(reshape)
library(e1071)

# DECIDE RUN MODE
DEBUG_MODE = FALSE

# SETUP OUTPUT DIRECTORY
if (!file.exists("outputs_statistical"))
{
  dir.create("outputs_statistical")
}

# LOAD HELPERS
source('freq_weight.R')
source('skill_statistics.R')
source('log_seq.R')
source('quantile_trim.R')

# KEY PARAMETERS
if (DEBUG_MODE==TRUE)
{
  CORES = 1
  REPLICATES = 1
  GRID_POINTS = 2
  MIN_POINTS = 1e1
  MAX_POINTS = 1e3
} else
{
  CORES = 16 # number of cores to parallel process on
  REPLICATES = 5
  GRID_POINTS = 20
  MIN_POINTS = 1e1
  MAX_POINTS = 1e4
}


###### MAIN FUNCTIONS
predict_rf <- function(yvar, assemblages, rows_train, method, num_species)
{
  print(yvar) # DEBUG
  
  if (length(rows_train)==0) # if there is no training data at the expected richness
  {
    return(NULL)
  }
  
  rows_test = setdiff(1:nrow(assemblages), rows_train)
  
  # select training subset
  data_for_rf_training = assemblages[rows_train,]

  if (method=='e2e + reshuffle')
  {
    # shuffle the x data
    data_for_rf_training[,names(data_for_rf_training)[1:num_species]] = data_for_rf_training[sample(1:nrow(data_for_rf_training)),names(data_for_rf_training)[1:num_species]]
  }

  if (yvar %in% c("_composition","_abundance"))
  {
    yvar_this = names(data_for_rf_training)[grep("star",names(data_for_rf_training))]
    
    if (yvar=="_composition")
    {
      # convert to presence/absence data
      data_for_rf_training[,yvar_this] = data_for_rf_training[,yvar_this] > 0
      for (var in yvar_this)
      {
        data_for_rf_training[,var] = factor(data_for_rf_training[,var],levels=c(FALSE,TRUE))
      }
    }
    else # if we have abundance
    {
      # convert the abundances to 
      ### BB 
      numeric_vals = data_for_rf_training[,yvar_this] %>% as.matrix %>% as.numeric
      quantiles = quantile(numeric_vals[numeric_vals>0], seq(0,1,length.out=9), na.rm=TRUE) # assume ten total classes
      maxval = max(numeric_vals, na.rm=TRUE)
      if (maxval==0)
      {
        maxval=1e-16 # a hack to get the breaks to work below
      }
      
      data_for_rf_training[,yvar_this] = data_for_rf_training[,yvar_this] %>% 
        mutate(across(everything(), function(x) {
          breaks = na.omit(as.numeric(unique(c(0, quantiles, maxval))))
          bin_means = (head(breaks,-1) + tail(breaks,-1))/2
          bin_means[1] = 0
          #print(bin_means) ### DEBUG
          return(cut(x, breaks=breaks,right=TRUE,include.lowest=TRUE,labels=bin_means))
        }
        ))
    }
    
    if (method %in% c('e2e','e2e + reshuffle'))
    {
      formula_this_multivar = formula(sprintf("Multivar(%s)~%s",paste(yvar_this,collapse=", "),paste(names(data_for_rf_training)[1:num_species],collapse="+")))
      
      m_rf_multivar = rfsrc(formula=formula_this_multivar,
                            data = data_for_rf_training,
                            forest = TRUE,
                            importance="none")
      
      values_predicted_raw = predict(object=m_rf_multivar, 
                                     newdata=assemblages[rows_test,1:num_species])
      if (yvar=="_composition")
      {
        values_predicted = as.data.frame(sapply(values_predicted_raw$classOutput, function(x) {as.logical(x$class)}, simplify=FALSE))
      }
      else # if abundance
      {
        # convert the class predictions back to numeric values
        values_predicted = as.data.frame(sapply(values_predicted_raw$classOutput, function(x) {x$class}, simplify=FALSE))
        values_predicted = values_predicted %>% mutate(across(everything(), function(x) {as.numeric(as.character(x))}))
      }
    }
    else if (method=='naive')
    {
      # use the input presence/absences
      if (yvar=="_composition")
      {
        values_predicted = assemblages[rows_test,1:num_species] > 0
      }
      else
      {
        values_predicted = assemblages[rows_test,1:num_species]
      }
      
      m_rf_multivar = NULL
    }
    
    if (yvar=="_composition")
    {
      values_observed = assemblages[rows_test,yvar_this] > 0
    }
    else # if abundance
    {
      values_observed = assemblages[rows_test,yvar_this]
    }
    
    # if we have composition
    if (yvar=="_composition")
    {
      balanced_accuracy = NA # do a casewise analysis
      try(balanced_accuracy <- balanced_accuracy_casewise_mean(pred=values_predicted, obs=values_observed))
      
      return(list(model=m_rf_multivar,
                  pred=values_predicted,
                  obs=values_observed,
                  ba=balanced_accuracy))
    }
    else # if we have abundance
    {
      mae_mean = mean_absolute_error_casewise_mean(pred=values_predicted, obs=values_observed)
      
      #plot(as.numeric(as.matrix(values_predicted)), as.numeric(as.matrix(values_observed))) # DEBUG
      #abline(0,1,col='red') # DEBUG
      
      final_result = list(model=m_rf_multivar,
                          pred=values_predicted,
                          obs=values_observed,
                          mae_mean=mae_mean)
      
      return(final_result)      
    }
  }
  else # if we are doing single variable prediction
  {
    if (method %in% c('e2e','e2e + reshuffle'))
    {
      # add weights
      formula_this_1var = formula(sprintf("%s~%s",yvar,paste(names(data_for_rf_training)[1:num_species],collapse="+")))
      
      case_weights =  freq_weight(as.numeric(data_for_rf_training[,yvar]))
      
      m_rf_1var = ranger(data=data_for_rf_training,
                         formula=formula_this_1var,
                         importance='permutation',
                         case.weights = case_weights,
                         verbose=TRUE)
      
      values_predicted = predict(m_rf_1var, data=assemblages[rows_test,])$predictions
      
    }
    else if (method=='naive') # naive approach
    {
      if (is.factor(data_for_rf_training[,yvar]))
      {
        # randomly sample values of the factor from the training data
        values_predicted = sample(x = data_for_rf_training[,yvar],size = length(rows_test),replace = TRUE)
      }
      else
      {
        if (yvar=="richness") # assume all species coexist
        {
          values_predicted = apply(assemblages[rows_test,1:num_species],1,sum) # use training data richness
        }
        else         # pick mean value of the continuous variable (e.g. richness)
        {
          values_predicted = mean(data_for_rf_training[,yvar])
        }
      }
      
      m_rf_1var = NULL
    }
    
    values_observed = assemblages[rows_test,yvar]
    
    results = data.frame(pred=values_predicted, obs=values_observed)
    
    if(is.factor(assemblages[,yvar]) | is.logical(assemblages[,yvar]))
    {
      confusion = confusionMatrix(factor(results$pred,levels=c(FALSE,TRUE)), factor(results$obs,levels=c(FALSE,TRUE)))
      
      balanced_accuracy = confusion$byClass["Balanced Accuracy"]
      
      return(list(model=m_rf_1var,
                  pred=values_predicted,
                  obs=values_observed,
                  cm=confusion,
                  ba=balanced_accuracy))
    }
    else
    {
      mae = mean_absolute_error(pred=values_predicted, obs=values_observed)

      return(list(model=m_rf_1var,
                  pred=values_predicted,
                  obs=values_observed,
                  mae=mae))
    }
  }
}




do_predictions <- function(input_file,fn,
                           num_species, 
                           num_replicates_in_data = 1, 
                           num_replicates_in_rf=REPLICATES,
                           num_grid_points=GRID_POINTS,
                           min_points=MIN_POINTS,
                           max_points=MAX_POINTS)
{
  data = input_file %>%
    mutate(feasible.and.stable = factor(stable & feasible)) %>%
    select(-feasible, -stable)
  
  # remove outliers for all datasets
  data = quantile_trim(data)
  
  # remove missing cases that arose from the above
  which_rows_na = data %>% 
    select(contains("star")) %>% 
    rowSums %>%
    is.na %>%
    which
  print(sprintf("Removed %d NA case rows",length(which_rows_na)))
  data = data[-which_rows_na,]

  # make the sample size sequence  
  sample_size_seq_all = unique(round(log_seq(min_points, max_points, length.out=num_grid_points)))
  # trim to only the sizes that are compatible with the dataset
  sample_size_seq_all = sample_size_seq_all[sample_size_seq_all <= nrow(data)]
  
  results_table = expand.grid(rep=1:num_replicates_in_rf, 
                              method=c('naive','e2e + reshuffle','e2e'),
                              richness_mae=NA,
                              feasible_and_stable_balanced_accuracy=NA,
                              composition_balanced_accuracy_mean=NA,
                              abundance_mae_mean=NA,
                              sampling_strategy=c("high-1","high-2","high-3",
                                                  "low-1","low-2","low-3",
                                                  "mixed"),
                              num_train=sample_size_seq_all,
                              num_species_dataset=NA,
                              num_replicates_dataset=NA,
                              num_cases_dataset=NA,
                              num_losses_mean_train=NA,
                              abundance_q95_dataset=NA,
                              abundance_skewness_dataset=NA,
                              abundance_skewness_nonzero_dataset=NA,
                              abundance_final_skewness_mean_train=NA,
                              abundance_final_skewness_nonzero_mean_train=NA
  )
  
  results_list = mclapply(1:nrow(results_table), function(i) # DEBUG #lapply(1:nrow(results_table), function(i)#
  {
    #cat('.')
    print(data.frame(file=fn, i=i, fraction_finished=i/nrow(results_table), results_table[i,])) #DEBUG

    # set sample size
    n_train = results_table$num_train[i]
    
    if (results_table$sampling_strategy[i]=="mixed")
    {
      if (n_train == nrow(data))
      {
        print('all points in training, no points, skipping')
        return(NULL)
      }
      else
      {
        rows_train = sample(x=1:nrow(data), size=n_train)
      }
    }
    else if (results_table$sampling_strategy[i] %in% c("low-1","low-2","low-3"))
    {
      max_richness = as.numeric(gsub("low-","",as.character(results_table$sampling_strategy[i]),fixed=TRUE))
      initial_richness = apply(data[,1:num_species],1,sum)
      
      which_rows = which(initial_richness <= max_richness)
      
      if (n_train > length(which_rows)) # if we are requesting more points than exist in the number of cases available
      {
        print('too many points requested for low-* sampling, skipping')
        return(NULL)
      }
      
      rows_train = sample(x=which_rows, size=n_train)
    }
    else if (results_table$sampling_strategy[i] %in% c("high-1","high-2","high-3"))
    {
      min_richness = num_species - as.numeric(gsub("high-","",as.character(results_table$sampling_strategy[i]),fixed=TRUE))
      initial_richness = apply(data[,1:num_species],1,sum)
      
      which_rows = which(initial_richness >= min_richness)
      
      if (n_train > length(which_rows)) # if we are requesting more points than exist in the number of cases available
      {
        print('too many points requested for high-* sampling, skipping')
        return(NULL)
      }
      
      rows_train = sample(x=which_rows, size=n_train)
    }
    
    prediction_abundance = predict_rf(yvar = "_abundance",
                                      assemblages = data,
                                      rows_train = rows_train,
                                      method = results_table$method[i],
                                      num_species=num_species)

    if(!is.null(prediction_abundance))
    {
      results_table$abundance_mae_mean[i]=prediction_abundance$mae_mean
      write.csv(prediction_abundance$pred, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_num_train=%d_sampling_strategy=%s_abundance_pred.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$sampling_strategy[i]), row.names=FALSE)
      write.csv(prediction_abundance$obs, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_num_train=%d_sampling_strategy=%s_abundance_obs.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$sampling_strategy[i]), row.names=FALSE)
    }   
    #print(results_table[i,,drop=FALSE]) # DEBUG
    
    prediction_composition = predict_rf(yvar = "_composition",
                                        assemblages = data,
                                        rows_train = rows_train,
                                        method = results_table$method[i],
                                        num_species=num_species)
    if (!is.null(prediction_composition))
    {
      results_table$composition_balanced_accuracy_mean[i]=prediction_composition$ba
      write.csv(prediction_composition$pred, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_num_train=%d_sampling_strategy=%s_composition_pred.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$sampling_strategy[i]), row.names=FALSE)
      write.csv(prediction_composition$obs, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_num_train=%d_sampling_strategy=%s_composition_obs.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$sampling_strategy[i]), row.names=FALSE)
    }
    #print(results_table[i,,drop=FALSE]) # DEBUG
    
    prediction_richness = predict_rf(yvar = 'richness',
                                     assemblages = data,
                                     rows_train = rows_train,
                                     method = results_table$method[i],
                                     num_species=num_species)
    if (!is.null(prediction_richness))
    {
      results_table$richness_mae[i]=prediction_richness$mae
      write.csv(prediction_richness$pred, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_num_train=%d_sampling_strategy=%s_richness_pred.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$sampling_strategy[i]), row.names=FALSE)
      write.csv(prediction_richness$obs, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_num_train=%d_sampling_strategy=%s_richness_obs.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$sampling_strategy[i]), row.names=FALSE)
    }
    #print(results_table[i,,drop=FALSE]) # DEBUG
    
    if(!all(is.na(data$feasible.and.stable)))
    {
      prediction_fs = predict_rf(yvar = 'feasible.and.stable',
                                 assemblages = data,
                                 rows_train = rows_train,
                                 method = results_table$method[i],
                                 num_species=num_species)
    }
    else
    {
      prediction_fs = NULL
    }
    if (!is.null(prediction_fs))
    {
      results_table$feasible_and_stable_balanced_accuracy[i]=prediction_fs$ba
      write.csv(prediction_fs$pred, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_num_train=%d_sampling_strategy=%s_fs_pred.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$sampling_strategy[i]), row.names=FALSE)
      write.csv(prediction_fs$obs, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_num_train=%d_sampling_strategy=%s_fs_obs.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$sampling_strategy[i]), row.names=FALSE)
    }
    #print(results_table[i,,drop=FALSE]) # DEBUG
    
    results_table$num_train[i] = n_train
    results_table$num_species_dataset[i] = num_species
    results_table$num_replicates_dataset[i] = num_replicates_in_data
    results_table$num_cases_dataset[i] = nrow(data)
    
    # get training set outcomes
    final_abundances_train = data[rows_train,] %>% select(contains("star")) %>% as.matrix
    initial_conditions_train = data[rows_train,1:num_species] %>% as.matrix
    
    # count # of species that were present but went absent
    results_table$num_losses_mean_train[i] = mean(apply((final_abundances_train==0) & (initial_conditions_train==1), 1, sum, na.rm=TRUE))

    # figure out abundance distribution in training
    results_table$abundance_final_skewness_mean_train[i] = skewness(as.numeric(final_abundances_train), na.rm=TRUE)
    results_table$abundance_final_skewness_nonzero_mean_train[i] = skewness(as.numeric(final_abundances_train)[as.numeric(final_abundances_train)>0], na.rm=TRUE)
    
    # determine the overall dataset 95% abundance quantile
    abundances_dataset_all = data %>% select(contains("star")) %>% as.matrix %>% as.numeric
    results_table$abundance_q95_dataset[i] = quantile(abundances_dataset_all, 0.95, na.rm=TRUE)

    # determine overall dataset skewness
    results_table$abundance_skewness_dataset[i] = skewness(abundances_dataset_all, na.rm=TRUE)
    results_table$abundance_skewness_nonzero_dataset[i] = skewness(abundances_dataset_all[abundances_dataset_all>0], na.rm=TRUE)

    return(results_table[i,,drop=FALSE])
  }, mc.cores=CORES)
  
  results_df = rbindlist(results_list)
  
  # write results table
  write.csv(results_df, file=sprintf('outputs_statistical/results_%s.csv',fn), row.names=FALSE)
}











########################
# DO ANALYSES
data_assemblages_M_11 = read.csv('data_glv/assemblages_M_11.csv')
do_predictions(data_assemblages_M_11,
               fn = 'mouse_gut',
               num_species = 11,
               num_replicates_in_data = 1)

data_soil_bacteria_8 = read.csv('data_friedman_gore/data_friedman_gore.csv')
do_predictions(data_soil_bacteria_8,
               fn = 'soil_bacteria',
               num_species = 8, 
               num_replicates_in_data = 2) # this is an underestimate but should not cause problems

data_assemblages_H_12 = read.csv('data_glv/assemblages_H_12.csv')
do_predictions(data_assemblages_H_12,
               'human_gut',
               num_species = 12,
               num_replicates_in_data = 1)

data_assemblages_glv_16 = read.csv('data_glv/assemblages_glv_16.csv')
data_assemblages_glv_16 = data_assemblages_glv_16 %>% sample_n(2^14)
do_predictions(data_assemblages_glv_16,
               fn = 'glv_simulated',
               num_species = 16,
               num_replicates_in_data = 1)

data_annual_plant_18 = read.csv('data_annual_plant/assemblages_annual_plant_18.csv') %>%
  mutate(stable=replace(stable, is.na(stable), 1)) # as there is one missing case (the all zeros case)
data_annual_plant_18 = data_annual_plant_18 %>% sample_n(2^14)
do_predictions(data_annual_plant_18,
               fn = 'annual_plant',
               num_species = 18,
               num_replicates_in_data = 1)

data_sortie_9_3 = read.csv('data_sortie/data_sortie.csv')
do_predictions(data_sortie_9_3,
               fn = 'sortie-nd_plants',
               num_species = 9,
               num_replicates_in_data = 3)

data_fly_5 = read.csv('data_fly/data_fly.csv')
do_predictions(data_fly_5,
               fn = 'fly_gut',
               num_species = 5,
               num_replicates_in_data = 48)

data_assemblages_cedar_creek_18 = read.csv('data_cedar_creek/cedar_creek_2018.csv')
do_predictions(data_assemblages_cedar_creek_18,
               fn = 'cedar_creek_plants',
               num_species = 18,
               num_replicates_in_data = 1)
