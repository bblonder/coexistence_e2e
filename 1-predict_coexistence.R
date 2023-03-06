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

# DECIDE RUN MODE
DEBUG_MODE = FALSE

# SETUP OUTPUT DIRECTORY
if (!file.exists("outputs_statistical"))
{
  dir.create("outputs_statistical")
}

# LOAD HELPERS
source('utils/freq_weight.R')
source('utils/skill_statistics.R')
source('utils/log_seq.R')
source('utils/quantile_trim.R')

# KEY PARAMETERS
if (DEBUG_MODE==TRUE)
{
  CORES = 1
  REPLICATES = 1
  GRID_POINTS = 1
  MIN_POINTS = 1e2
  MAX_POINTS = 1e2
} else
{
  CORES = 4
  REPLICATES = 10
  GRID_POINTS = 20
  MIN_POINTS = 1e1
  MAX_POINTS = 1e4
}

###### MAIN FUNCTIONS
predict_model <- function(yvar, assemblages, rows_train, method, num_species)
{
  cat(yvar) # DEBUG
  cat(' ')
  
  if (length(rows_train)==0) # if there is no training data at the expected richness
  {
    return(NULL)
  }
  
  rows_test = setdiff(1:nrow(assemblages), rows_train)
  
  # select training subset
  data_for_rf_training = assemblages[rows_train,]

  if (method=='love + reshuffle')
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
    
    if (method %in% c('love','love + reshuffle'))
    {
      formula_this_multivar = formula(sprintf("Multivar(%s)~%s",paste(yvar_this,collapse=", "),paste(names(data_for_rf_training)[1:num_species],collapse="+")))
      
      m_rf_multivar = rfsrc(formula=formula_this_multivar,
                            data = data_for_rf_training,
                            forest = TRUE,
                            importance="none",
                            num.trees=500,
                            mtry=ceiling(sqrt(num_species)),
                            min.node.size=ceiling(sqrt(num_species))
                            )
      
      values_predicted_raw_test = predict(object=m_rf_multivar, 
                                     newdata=assemblages[rows_test,1:num_species])
      values_predicted_raw_train = predict(object=m_rf_multivar, 
                                          newdata=assemblages[rows_train,1:num_species])
      if (yvar=="_composition")
      {
        values_predicted_test = as.data.frame(sapply(values_predicted_raw_test$classOutput, function(x) {as.logical(x$class)}, simplify=FALSE))
        values_predicted_train = as.data.frame(sapply(values_predicted_raw_train$classOutput, function(x) {as.logical(x$class)}, simplify=FALSE))
      }
      else # if abundance
      {
        # convert the class predictions back to numeric values
        values_predicted_test = as.data.frame(sapply(values_predicted_raw_test$classOutput, function(x) {x$class}, simplify=FALSE))
        values_predicted_test = values_predicted_test %>% mutate(across(everything(), function(x) {as.numeric(as.character(x))}))
        
        values_predicted_train = as.data.frame(sapply(values_predicted_raw_train$classOutput, function(x) {x$class}, simplify=FALSE))
        values_predicted_train = values_predicted_train %>% mutate(across(everything(), function(x) {as.numeric(as.character(x))}))
      }
    }
    else if (method=='naive')
    {
      # if composition use the input presence/absences
      if (yvar=="_composition")
      {
        values_predicted_test = assemblages[rows_test,1:num_species] > 0
        values_predicted_train = assemblages[rows_train,1:num_species] > 0
      }
      # if abundance use the mean training values masked by the input presence/absences
      else
      {
        values_predicted_test = assemblages[rows_test,1:num_species] * colMeans(assemblages[rows_train,yvar_this])
        values_predicted_train = assemblages[rows_train,1:num_species] * colMeans(assemblages[rows_train,yvar_this])
      }
      
      m_rf_multivar = NULL
    }
    
    if (yvar=="_composition")
    {
      values_observed_test = assemblages[rows_test,yvar_this] > 0
      values_observed_train = assemblages[rows_train,yvar_this] > 0
    }
    else # if abundance
    {
      values_observed_test = assemblages[rows_test,yvar_this]
      values_observed_train = assemblages[rows_train,yvar_this]
    }
    
    # if we have composition
    if (yvar=="_composition")
    {
      balanced_accuracy_test = NA # do a casewise analysis
      try(balanced_accuracy_test <- balanced_accuracy_casewise_mean(pred=values_predicted_test, obs=values_observed_test))
      
      balanced_accuracy_train = NA # do a casewise analysis
      try(balanced_accuracy_train <- balanced_accuracy_casewise_mean(pred=values_predicted_train, obs=values_observed_train))
      
      return(list(model=m_rf_multivar,
                  pred_test=values_predicted_test,
                  obs_test=values_observed_test,
                  ba_test=balanced_accuracy_test,
                  pred_train=values_predicted_train,
                  obs_train=values_observed_train,
                  ba_train=balanced_accuracy_train))
    }
    else # if we have abundance
    {
      mae_mean_test = mean_absolute_error_casewise_mean(pred=values_predicted_test, obs=values_observed_test)
      mae_mean_train = mean_absolute_error_casewise_mean(pred=values_predicted_train, obs=values_observed_train)
      
      #plot(as.numeric(as.matrix(values_predicted)), as.numeric(as.matrix(values_observed))) # DEBUG
      #abline(0,1,col='red') # DEBUG
      
      final_result = list(model=m_rf_multivar,
                          pred_test=values_predicted_test,
                          obs_test=values_observed_test,
                          mae_mean_test=mae_mean_test,
                          pred_train=values_predicted_train,
                          obs_train=values_observed_train,
                          mae_mean_train=mae_mean_train)
      
      return(final_result)      
    }
  }
  else # if we are doing single variable prediction
  {
    if (method %in% c('love','love + reshuffle'))
    {
      # add weights
      formula_this_1var = formula(sprintf("%s~%s",yvar,paste(names(data_for_rf_training)[1:num_species],collapse="+")))
      
      case_weights =  freq_weight(as.numeric(data_for_rf_training[,yvar]))
      
      m_rf_1var = ranger(data=data_for_rf_training,
                         formula=formula_this_1var,
                         importance='permutation',
                         case.weights = case_weights,
                         verbose=TRUE,
                         num.trees=500,
                         mtry=ceiling(sqrt(num_species)),
                         min.node.size=ceiling(sqrt(num_species))
                         )
      
      values_predicted_test = predict(m_rf_1var, data=assemblages[rows_test,])$predictions
      values_predicted_train = predict(m_rf_1var, data=assemblages[rows_train,])$predictions
    }
    else if (method=='naive') # naive approach
    {
      if (is.factor(data_for_rf_training[,yvar]))
      {
        # randomly sample values of the factor from the training data
        values_predicted_test = sample(x = data_for_rf_training[,yvar],size = length(rows_test),replace = TRUE)
        values_predicted_train = sample(x = data_for_rf_training[,yvar],size = length(rows_train),replace = TRUE)
      }
      else
      {
        if (yvar=="richness") # assume all species coexist
        {
          values_predicted_test = apply(assemblages[rows_test,1:num_species],1,sum) # use training data richness
          values_predicted_train = apply(assemblages[rows_train,1:num_species],1,sum) # use training data richness
        }
        else         # pick mean value of the continuous variable (e.g. richness)
        {
          values_predicted_test = mean(data_for_rf_training[,yvar])
          values_predicted_train = mean(data_for_rf_training[,yvar])
        }
      }
      
      m_rf_1var = NULL
    }
    
    values_observed_test = assemblages[rows_test,yvar]
    values_observed_train = assemblages[rows_train,yvar]
    
    results_test = data.frame(pred=values_predicted_test, obs=values_observed_test)
    results_train = data.frame(pred=values_predicted_train, obs=values_observed_train)
    
    if(is.factor(assemblages[,yvar]) | is.logical(assemblages[,yvar]))
    {
      confusion_test = confusionMatrix(factor(results_test$pred,levels=c(FALSE,TRUE)), factor(results_test$obs,levels=c(FALSE,TRUE)))
      confusion_train = confusionMatrix(factor(results_train$pred,levels=c(FALSE,TRUE)), factor(results_train$obs,levels=c(FALSE,TRUE)))
      
      balanced_accuracy_test = confusion_test$byClass["Balanced Accuracy"]
      balanced_accuracy_train = confusion_train$byClass["Balanced Accuracy"]
      
      return(list(model=m_rf_1var,
                  pred_test=values_predicted_test,
                  obs_test=values_observed_test,
                  cm_test=confusion_test,
                  ba_test=balanced_accuracy_test,
                  pred_train=values_predicted_train,
                  obs_train=values_observed_train,
                  cm_train=confusion_train,
                  ba_train=balanced_accuracy_train))
    }
    else
    {
      mae_test = mean_absolute_error(pred=values_predicted_test, obs=values_observed_test)
      mae_train = mean_absolute_error(pred=values_predicted_train, obs=values_observed_train)

      return(list(model=m_rf_1var,
                  pred_test=values_predicted_test,
                  obs_test=values_observed_test,
                  mae_test=mae_test,
                  pred_train=values_predicted_train,
                  obs_train=values_observed_train,
                  mae_train=mae_train))
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
    mutate(feasible.and.stable = factor(stable & feasible, levels=c(FALSE, TRUE))) %>%
    select(-feasible, -stable)
  
  # flag quantile outliers
  data = quantile_max_trim(data)
  
  # remove missing cases that arose from the above
  which_rows_na = data %>% 
    select(contains("star")) %>% 
    rowSums %>%
    is.na %>%
    which
  print(sprintf("Removed %d problematic case rows",length(which_rows_na)))
  data = data[-which_rows_na,]

  # make the sample size sequence  
  sample_size_seq_all = unique(round(log_seq(min_points, max_points, length.out=num_grid_points)))
  # trim to only the sizes that are compatible with the dataset
  sample_size_seq_all = sample_size_seq_all[sample_size_seq_all <= nrow(data)]
  
  results_table = expand.grid(rep=1:num_replicates_in_rf, 
                              method=c('naive','love + reshuffle','love'),
                              richness_mae_test=NA,
                              feasible_and_stable_balanced_accuracy_test=NA,
                              composition_balanced_accuracy_mean_test=NA,
                              abundance_mae_mean_test=NA,
                              richness_mae_train=NA,
                              feasible_and_stable_balanced_accuracy_train=NA,
                              composition_balanced_accuracy_mean_train=NA,
                              abundance_mae_mean_train=NA,
                              experimental_design=c("high-1","high-2",
                                                  "low-2","low-3",
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
  
  indices = 1:nrow(results_table)
  
  results_list = mclapply(indices, function(i)
  {
    #cat('.')
    cat(sprintf('%d %d', i, nrow(results_table)))
    cat('\n')
    #print(data.frame(file=fn, i=i, fraction_finished=i/nrow(results_table), results_table[i,])) #DEBUG

    # set sample size
    n_train = results_table$num_train[i]
    
    if (results_table$experimental_design[i]=="mixed")
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
    else if (results_table$experimental_design[i] %in% c("low-1","low-2","low-3"))
    {
      max_richness = as.numeric(gsub("low-","",as.character(results_table$experimental_design[i]),fixed=TRUE))
      initial_richness = apply(data[,1:num_species],1,sum)
      
      which_rows = which(initial_richness <= max_richness)
      
      if (n_train > length(which_rows)) # if we are requesting more points than exist in the number of cases available
      {
        print('too many points requested for low-* sampling, skipping')
        return(NULL)
      }
      
      rows_train = sample(x=which_rows, size=n_train)
    }
    else if (results_table$experimental_design[i] %in% c("high-1","high-2","high-3"))
    {
      min_richness = num_species - as.numeric(gsub("high-","",as.character(results_table$experimental_design[i]),fixed=TRUE))
      initial_richness = apply(data[,1:num_species],1,sum)
      
      which_rows = which(initial_richness >= min_richness)
      
      if (n_train > length(which_rows)) # if we are requesting more points than exist in the number of cases available
      {
        print('too many points requested for high-* sampling, skipping')
        return(NULL)
      }
      
      rows_train = sample(x=which_rows, size=n_train)
    }
    
    # write out the experimental conditions if there are enough training rows
    try(write.csv(data[rows_train,], 
              file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=abundance_output=experiment_train.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE))
    try(write.csv(data[setdiff(1:nrow(data),rows_train),], 
              file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=abundance_output=experiment_test.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE))
        
    prediction_abundance = predict_model(yvar = "_abundance",
                                      assemblages = data,
                                      rows_train = rows_train,
                                      method = results_table$method[i],
                                      num_species=num_species)

    if(!is.null(prediction_abundance))
    {
      results_table$abundance_mae_mean_test[i]=prediction_abundance$mae_mean_test
      results_table$abundance_mae_mean_train[i]=prediction_abundance$mae_mean_train
      write.csv(prediction_abundance$pred_test, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=abundance_output=pred_test.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
      write.csv(prediction_abundance$obs_test, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=abundance_output=obs_test.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
      write.csv(prediction_abundance$pred_train, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=abundance_output=pred_train.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
      write.csv(prediction_abundance$obs_train, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=abundance_output=obs_train.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
      
    }   
    #print(results_table[i,,drop=FALSE]) # DEBUG
    
    prediction_composition = predict_model(yvar = "_composition",
                                        assemblages = data,
                                        rows_train = rows_train,
                                        method = results_table$method[i],
                                        num_species=num_species)
    if (!is.null(prediction_composition))
    {
      results_table$composition_balanced_accuracy_mean_test[i]=prediction_composition$ba_test
      results_table$composition_balanced_accuracy_mean_train[i]=prediction_composition$ba_train
      write.csv(prediction_composition$pred_test, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=composition_output=pred_test.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
      write.csv(prediction_composition$obs_test, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=composition_output=obs_test.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
      write.csv(prediction_composition$pred_train, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=composition_output=pred_train.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
      write.csv(prediction_composition$obs_train, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=composition_output=obs_train.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
    }
    #print(results_table[i,,drop=FALSE]) # DEBUG
    
    prediction_richness = predict_model(yvar = 'richness',
                                     assemblages = data,
                                     rows_train = rows_train,
                                     method = results_table$method[i],
                                     num_species=num_species)
    if (!is.null(prediction_richness))
    {
      results_table$richness_mae_test[i]=prediction_richness$mae_test
      results_table$richness_mae_train[i]=prediction_richness$mae_train
      write.csv(prediction_richness$pred_test, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=richness_output=pred_test.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
      write.csv(prediction_richness$obs_test, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=richness_output=obs_test.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
      write.csv(prediction_richness$pred_train, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=richness_output=pred_train.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
      write.csv(prediction_richness$obs_train, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=richness_output=obs_train.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
    }
    #print(results_table[i,,drop=FALSE]) # DEBUG
    
    # if(!all(is.na(data$feasible.and.stable)))
    # {
    #   prediction_fs = predict_model(yvar = 'feasible.and.stable',
    #                              assemblages = data,
    #                              rows_train = rows_train,
    #                              method = results_table$method[i],
    #                              num_species=num_species)
    # }
    # else
    # {
    #   prediction_fs = NULL
    # }
    # if (!is.null(prediction_fs))
    # {
    #   results_table$feasible_and_stable_balanced_accuracy[i]=prediction_fs$ba
    #   write.csv(prediction_fs$pred, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=fs_output=pred.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
    #   write.csv(prediction_fs$obs, file=sprintf('outputs/statistical/table_dataset=%s_i=%d_method=%s_rep=%d_num_train=%d_experimental_design=%s_response_var=fs_output=obs.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$num_train[i], results_table$experimental_design[i]), row.names=FALSE)
    # }
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
    
    cat('\n')
    
    return(results_table[i,,drop=FALSE])
  }, mc.cores=CORES)
  
  # write results table
  results_df = NULL
  try(results_df <- rbindlist(results_list))
  if (!is.null(results_df))
  {
    write.csv(results_df, file=sprintf('outputs/statistical/results_%s.csv',fn), row.names=FALSE)
  }
  # save the raw output too in case of a rbind issue for error checking
  saveRDS(results_list,file=sprintf('outputs/statistical/results_%s.Rdata',fn))
  
  indices_errors = which(sapply(results_list,class)=="try-error")
  if (length(indices_errors) > 0)
  {
    print(data.frame(indices_errors, sapply(indices_errors, function(j) { attr(results_list[[j]],"condition") })))
  }
  
  return(results_list)  
}











########################
# DO ANALYSES
set.seed(1)
data_annual_plant_18 = read.csv('data/annual_plant/assemblages_annual_plant_18.csv')
#if (DEBUG_MODE==TRUE)
#{
#  data_annual_plant_18 = data_annual_plant_18 %>% sample_n(2^14)
#}
do_predictions(data_annual_plant_18,
               fn = 'annual_plant',
               num_species = 18,
               num_replicates_in_data = 1)

set.seed(1)
data_assemblages_M_11 = read.csv('data/glv/assemblages_M_11.csv')
do_predictions(data_assemblages_M_11,
               fn = 'mouse_gut',
               num_species = 11,
               num_replicates_in_data = 1)

set.seed(1)
data_soil_bacteria_8 = read.csv('data/friedman_gore/data_friedman_gore.csv')
do_predictions(data_soil_bacteria_8,
               fn = 'soil_bacteria',
               num_species = 8, 
               num_replicates_in_data = 2) # this is an underestimate but should not cause problems

set.seed(1)
data_assemblages_H_12 = read.csv('data/glv/assemblages_H_12.csv')
do_predictions(data_assemblages_H_12,
               'human_gut',
               num_species = 12,
               num_replicates_in_data = 1)

set.seed(1)
# data_assemblages_glv_16 = read.csv('data/glv/assemblages_glv_16.csv')
# #if (DEBUG_MODE==TRUE)
# #{
#   data_assemblages_glv_16 = data_assemblages_glv_16 %>% sample_n(2^14)
# #}
# do_predictions(data_assemblages_glv_16,
#                fn = 'glv_simulated',
#                num_species = 16,
#                num_replicates_in_data = 1)

set.seed(1)
data_sortie_9_3 = read.csv('data/sortie/data_sortie.csv')
do_predictions(data_sortie_9_3,
               fn = 'sortie-nd_plants',
               num_species = 9,
               num_replicates_in_data = 3)

set.seed(1)
data_fly_5 = read.csv('data/fly/data_fly.csv')
do_predictions(data_fly_5,
               fn = 'fly_gut',
               num_species = 5,
               num_replicates_in_data = 48)

set.seed(1)
data_assemblages_cedar_creek_18 = read.csv('data/cedar_creek/cedar_creek_2018.csv')
do_predictions(data_assemblages_cedar_creek_18,
               fn = 'cedar_creek_plants',
               num_species = 18,
               num_replicates_in_data = 1)
