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
source('casewise.R')
source('lseq.R')
source('quantile_trim.R')

# KEY PARAMETERS
if (DEBUG_MODE==TRUE)
{
  CORES = 1#16 # number of cores to parallel process on
  REPLICATES = 1#5
  GRID_POINTS = 1#11
  MIN_FRAC = 1e-2#1e-5
} else
{
  CORES = 16 # number of cores to parallel process on
  REPLICATES = 5
  GRID_POINTS = 11
  MIN_FRAC = 1e-5  
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
      quantiles = quantile(numeric_vals[numeric_vals>0], seq(0,1,length.out=(10-2)), na.rm=TRUE) # assume ten total classes
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
        values_predicted = as.data.frame(sapply(values_predicted_raw$classOutput, function(x) {as.logical(x$class)}))
      }
      else # if abundance
      {
        # convert the class predictions back to numeric values
        values_predicted = as.data.frame(sapply(values_predicted_raw$classOutput, function(x) {x$class}))
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
      r2 = NA # do a casewise analysis
      try(r2 <- r2_casewise_mean(pred=values_predicted, obs=values_observed))
      
      #plot(as.numeric(as.matrix(values_predicted)), as.numeric(as.matrix(values_observed))) # DEBUG
      #abline(0,1,col='red') # DEBUG
      
      final_result = list(model=m_rf_multivar,
                          pred=values_predicted,
                          obs=values_observed,
                          r2=r2)
      
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
    
    if(is.factor(assemblages[,yvar]))
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
      m_lm = NULL
      try(m_lm <- lm(pred~obs,data=results))
      if (!is.null(m_lm))
      {
        suppressWarnings(r2 <- summary(m_lm)$r.squared) # disable warnings for the perfect fits this may generate
      }
      else
      {
        r2 = NA
      }
      return(list(model=m_rf_1var,
                  pred=values_predicted,
                  obs=values_observed,
                  r2=r2))
    }
  }
}




do_predictions <- function(input_file,fn,
                           num_species, 
                           num_replicates_in_data = 1, 
                           num_replicates_in_rf=REPLICATES,
                           num_grid_points=GRID_POINTS,
                           min_frac = MIN_FRAC)
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
  
  # determine if we have all the possible cases
  dataset_complete = (num_replicates_in_data * 2^num_species == nrow(data))
  
  seq_all = lseq(min_frac,1,length.out=num_grid_points)
  
  if(dataset_complete==TRUE)
  {
    frac_this = seq_all
  }
  else
  {
    num_cases_this = nrow(input_file)
    frac_this = seq_all[which(seq_all <= num_cases_this/(num_replicates_in_data * 2^num_species) )]
  }
  
  results_table = expand.grid(rep=1:num_replicates_in_rf, 
                              method=c('naive','e2e + reshuffle','e2e'),
                              frac=frac_this,
                              richness.r2=NA,
                              feasible.and.stable.balanced_accuracy=NA,
                              composition.balanced_accuracy=NA,
                              abundance.r2=NA,
                              sampling_strategy=c("high-1","high-2","high-3",
                                                  "low-1","low-2","low-3",
                                                  "mixed"),
                              num_train=NA,
                              num_species=NA,
                              num_replicates_in_data=NA,
                              num_cases=NA,
                              richness_initial_mean=NA,
                              richness_final_mean=NA,
                              abundance_final_skewness_mean=NA,
                              abundance_final_skewness_nonzero_mean=NA
  )
  
  results_list = mclapply(1:nrow(results_table), function(i) # DEBUG #lapply(1:nrow(results_table), function(i)#
  {
    #cat('.')
    print(data.frame(file=fn, i=i, fraction_finished=i/nrow(results_table))) #DEBUG

    if (results_table$sampling_strategy[i]=="mixed")
    {
      n_train = ceiling(nrow(data)*results_table$frac[i])
      
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
      
      n_train = ceiling(nrow(data)*results_table$frac[i])
      if (n_train > length(which_rows)) # if we are requesting more points than exist in the number of cases available
      {
        print('too many points requested for low sampling, skipping')
        return(NULL)
      }
      
      rows_train = sample(x=which_rows, size=n_train)
    }
    else if (results_table$sampling_strategy[i] %in% c("high-1","high-2","high-3"))
    {
      min_richness = num_species - as.numeric(gsub("high-","",as.character(results_table$sampling_strategy[i]),fixed=TRUE))
      initial_richness = apply(data[,1:num_species],1,sum)
      
      which_rows = which(initial_richness >= min_richness)
      
      n_train = ceiling(nrow(data)*results_table$frac[i])
      if (n_train > length(which_rows)) # if we are requesting more points than exist in the number of cases available
      {
        print('too many points requested for high sampling, skipping')
        return(NULL)
      }
      
      rows_train = sample(x=which_rows, size=n_train)
    }
    #print(rows_train) # DEBUG
    
    prediction_abundance = predict_rf(yvar = "_abundance",
                                      assemblages = data,
                                      rows_train = rows_train,
                                      method = results_table$method[i],
                                      num_species=num_species)

    if(!is.null(prediction_abundance))
    {
      results_table$abundance.r2[i]=prediction_abundance$r2
      if (results_table$rep[i]==1) # just to cut down on duplicate outputs
      {
        write.csv(prediction_abundance$pred, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_frac=%f_sampling_strategy=%s_abundance_pred.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$frac[i], results_table$sampling_strategy[i]), row.names=FALSE)
        write.csv(prediction_abundance$obs, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_frac=%f_sampling_strategy=%s_abundance_obs.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$frac[i], results_table$sampling_strategy[i]), row.names=FALSE)
      }
    }   
    #print(results_table[i,,drop=FALSE]) # DEBUG
    
    prediction_composition = predict_rf(yvar = "_composition",
                                        assemblages = data,
                                        rows_train = rows_train,
                                        method = results_table$method[i],
                                        num_species=num_species)
    if (!is.null(prediction_composition))
    {
      results_table$composition.balanced_accuracy[i]=prediction_composition$ba
      write.csv(prediction_composition$pred, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_frac=%f_sampling_strategy=%s_composition_pred.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$frac[i], results_table$sampling_strategy[i]), row.names=FALSE)
      write.csv(prediction_composition$obs, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_frac=%f_sampling_strategy=%s_composition_obs.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$frac[i], results_table$sampling_strategy[i]), row.names=FALSE)
    }
    #print(results_table[i,,drop=FALSE]) # DEBUG
    
    prediction_richness = predict_rf(yvar = 'richness',
                                     assemblages = data,
                                     rows_train = rows_train,
                                     method = results_table$method[i],
                                     num_species=num_species)
    if (!is.null(prediction_richness))
    {
      results_table$richness.r2[i]=prediction_richness$r2
      write.csv(prediction_richness$pred, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_frac=%f_sampling_strategy=%s_richness_pred.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$frac[i], results_table$sampling_strategy[i]), row.names=FALSE)
      write.csv(prediction_richness$obs, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_frac=%f_sampling_strategy=%s_richness_obs.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$frac[i], results_table$sampling_strategy[i]), row.names=FALSE)
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
      results_table$feasible.and.stable.balanced_accuracy[i]=prediction_fs$ba
      write.csv(prediction_fs$pred, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_frac=%f_sampling_strategy=%s_fs_pred.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$frac[i], results_table$sampling_strategy[i]), row.names=FALSE)
      write.csv(prediction_fs$obs, file=sprintf('outputs_statistical/table_fn=%s_i=%d_method=%s_rep=%d_frac=%f_sampling_strategy=%s_fs_obs.csv', fn, i, results_table$method[i], results_table$rep[i], results_table$frac[i], results_table$sampling_strategy[i]), row.names=FALSE)
    }
    #print(results_table[i,,drop=FALSE]) # DEBUG
    
    
    results_table$num_train[i] = n_train
    
    results_table$num_species[i] = num_species
    results_table$num_replicates_in_data[i] = num_replicates_in_data
    results_table$num_cases[i] = nrow(data)
    
    results_table$richness_initial_mean[i] = mean(apply(data[rows_train,1:num_species],1,sum))
    results_table$richness_final_mean[i] = mean(apply(data[rows_train,names(data)[grep("star",names(data))] ] > 0,1,sum))
    values_abundance = as.numeric(as.matrix(prediction_abundance$obs))
    results_table$abundance_final_skewness_mean[i] = skewness(values_abundance, na.rm=TRUE)
    results_table$abundance_final_skewness_nonzero_mean[i] = skewness(values_abundance[values_abundance>0], na.rm=TRUE)
    

    #write.csv(results_table[i,,drop=FALSE],file=sprintf('temp_%s_%d.csv',fn,i),row.names=FALSE) # DEBUG
    
    return(results_table[i,,drop=FALSE])
  }, mc.cores=CORES) # DEBUG
  
  #saveRDS(results_list, file=sprintf('temp_results_%s.Rdata',fn)) #DEBUG
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

if (DEBUG_MODE==FALSE)
{
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
  
  data_soil_bacteria_8 = read.csv('data_friedman_gore/data_friedman_gore.csv')
  do_predictions(data_soil_bacteria_8,
                 fn = 'soil_bacteria',
                 num_species = 8, 
                 num_replicates_in_data = 2) # this is an underestimate but should not cause problems
  
  data_assemblages_cedar_creek_18 = read.csv('data_cedar_creek/cedar_creek_2018.csv')
  do_predictions(data_assemblages_cedar_creek_18,
                 fn = 'cedar_creek_plants',
                 num_species = 18,
                 num_replicates_in_data = 1)
}