library(ranger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(caret)
library(randomForestSRC)
library(parallel)
library(data.table)

freq_weight <- function(x)
{
  x_cut = cut(x, breaks=10)
  w = 1/table(x_cut)
  w = as.numeric(w)
  w[!is.finite(w)] = NA
  w <- w/sum(w,na.rm=T)
  names(w) = levels(x_cut)
  
  weights_by_x = as.numeric(w[x_cut])
  
  return(weights_by_x)
}

predict_rf <- function(yvar, assemblages, rows_train, method, num_species=log2(nrow(assemblages)))
{
  print(yvar) # DEBUG
  
  rows_test = setdiff(1:nrow(assemblages), rows_train)
  
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
    else
    {
      # otherwise leave as abundance
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
      else
      {
        values_predicted = as.data.frame(sapply(values_predicted_raw$regrOutput, function(x) {x$predicted}))
        
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
    else
    {
      values_observed = assemblages[rows_test,yvar_this]
    }
    
    if (yvar=="_composition")
    {
      confusion = confusionMatrix(factor(as.numeric(as.matrix(values_predicted)),levels=c(0,1)),factor(as.numeric(as.matrix(values_observed)),levels=c(0,1)))
      
      balanced_accuracy = confusion$byClass["Balanced Accuracy"]
      
      return(list(model=m_rf_multivar,
                  pred=values_predicted, 
                  obs=values_observed,
                  ba=balanced_accuracy))
    }
    else
    {
      results = data.frame(pred = as.numeric(as.matrix(values_predicted)), obs = as.numeric(as.matrix(values_observed)))
      message('across whole dataset not individual cases')
      m_lm = lm(pred~obs,data=results)

      intercept = coef(m_lm)[1]
      slope = coef(m_lm)[2]
      r2 = summary(m_lm)$r.squared
      
      final_result = list(model=m_rf_multivar,
                          df=results,
                          intercept=intercept,
                          slope=slope,
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
      
      values_predicted = as.numeric(predict(m_rf_1var, data=assemblages[rows_test,])$predictions)
      
    }
    else if (method=='naive') # naive approach
    {
      if (is.factor(data_for_rf_training[,yvar]))
      {
        # randomly sample values of the factor from the training data
        values_predicted = sample(x = as.numeric(data_for_rf_training[,yvar]),size = length(rows_test),replace = TRUE)
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
    
    values_observed = as.numeric(assemblages[rows_test,yvar])
    
    results = data.frame(pred=values_predicted, obs=values_observed)
    
    if(is.factor(assemblages[,yvar]))
    {
      confusion = confusionMatrix(factor(results$pred,levels=c(1,2)), factor(results$obs,levels=c(1,2)))
      
      balanced_accuracy = confusion$byClass["Balanced Accuracy"]
      
      return(list(model=m_rf_1var,
                  df=results,
                  cm=confusion,
                  ba=balanced_accuracy))
    }
    else
    {
      m_lm = lm(pred~obs,data=results)
      
      intercept = coef(m_lm)[1]
      slope = coef(m_lm)[2]
      r2 = summary(m_lm)$r.squared
      return(list(model=m_rf_1var,
                  df=results,
                  intercept=intercept,
                  slope=slope,
                  r2=r2))
    }
  }
}

do_predictions <- function(input_file,fn,
                           num_species, 
                           num_replicates_in_data = 1, 
                           dataset_complete, 
                           num_replicates_in_rf=5)
{
  data = input_file %>%
    mutate(feasible.and.stable = factor(stable & feasible))
  
  if(dataset_complete==TRUE)
  {
    frac_this = c(0.001,0.01,0.1,0.5,0.9)
  }
  else
  {
    num_cases_this = nrow(input_file)
    num_cases_total_this = 2^num_species * num_replicates_in_data
    frac_this = seq(1, num_cases_this, length.out=5) / num_cases_total_this
  }
  
  results_table = expand.grid(rep=1:num_replicates_in_rf, 
                              method=c('naive','e2e + reshuffle','e2e'),
                              frac=frac_this,
                              richness.r2=NA,
                              feasible.and.stable.balanced_accuracy=NA,
                              composition.balanced_accuracy=NA,
                              abundance.r2=NA
  )
  
  results_list = mclapply(1:nrow(results_table), function(i) # DEBUG #lapply(1:nrow(results_table), function(i)#
  {
    print(cbind(file=fn, i=i, frac=i/nrow(results_table), results_table[i,])) #DEBUG
    
    rows_train = sample(x=1:nrow(data), size=ceiling(nrow(data)*results_table$frac[i]))
    
    prediction_abundance = predict_rf(yvar = "_abundance",
                                      assemblages = data,
                                      rows_train = rows_train,
                                      method = results_table$method[i],
                                      num_species=num_species)
    results_table$abundance.r2[i]=prediction_abundance$r2
    print(results_table[i,,drop=FALSE]) # DEBUG
    
    prediction_composition = predict_rf(yvar = "_composition",
                                        assemblages = data,
                                        rows_train = rows_train,
                                        method = results_table$method[i],
                                        num_species=num_species)
    results_table$composition.balanced_accuracy[i]=prediction_composition$ba
    print(results_table[i,,drop=FALSE]) # DEBUG
    
    prediction_richness = predict_rf(yvar = 'richness',
                                     assemblages = data,
                                     rows_train = rows_train,
                                     method = results_table$method[i],
                                     num_species=num_species)
    results_table$richness.r2[i]=prediction_richness$r2
    print(results_table[i,,drop=FALSE]) # DEBUG
    
    prediction_fs = predict_rf(yvar = 'feasible.and.stable',
                               assemblages = data,
                               rows_train = rows_train,
                               method = results_table$method[i],
                               num_species=num_species)
    results_table$feasible.and.stable.balanced_accuracy[i]=prediction_fs$ba
    print(results_table[i,,drop=FALSE]) # DEBUG
    
    

    
    
    
    
    

    #write.csv(results_table[i,,drop=FALSE],file=sprintf('temp_%s_%d.csv',fn,i),row.names=FALSE) # DEBUG
    
    return(results_table[i,,drop=FALSE])
  }, mc.cores=2) # DEBUG
  
  #saveRDS(results_list, file=sprintf('results_%s.Rdata',fn)) #DEBUG
  results_df = rbindlist(results_list)
  
  
  # write results table
  write.csv(results_df, file=sprintf('outputs_statistical/results_%s.csv',fn), row.names=FALSE)
}












# do analyses
data_assemblages_cedar_creek_18 = read.csv('data_cedar_creek/cedar_creek_2018.csv')
do_predictions(data_assemblages_cedar_creek_18,
               fn = 'cedar_creek_n=18',
               num_species = 18,
               dataset_complete = FALSE)

data_sortie_9_3 = read.csv('data_sortie/data_sortie.csv')
do_predictions(data_sortie_9_3,
               fn = 'sortie_n=9_3',
               num_species = 9,
               dataset_complete = TRUE)

data_assemblages_H_12 = read.csv('data_glv/assemblages_H_12.csv')
do_predictions(data_assemblages_H_12,
               'glv_human_n=12',
               num_species = 12,
               dataset_complete = TRUE)

data_assemblages_M_11 = read.csv('data_glv/assemblages_M_11.csv')
do_predictions(data_assemblages_M_11,
               fn = 'glv_mouse_n=11',
               num_species = 11,
               dataset_complete = TRUE)

data_assemblages_glv_16 = read.csv('data_glv/assemblages_glv_16.csv')
data_assemblages_glv_16 = data_assemblages_glv_16 %>% sample_n(2^14)
do_predictions(data_assemblages_glv_16,
                fn = 'glv_simulated_n=16',
                num_species = 16, 
                dataset_complete = FALSE)

data_annual_plant_18 = read.csv('data_annual_plant/assemblages_annual_plant_18.csv') %>%
  mutate(stable=replace(stable, is.na(stable), 1)) # as there are missing cases
data_annual_plant_18 = data_annual_plant_18 %>% sample_n(2^14)
do_predictions(data_annual_plant_18,
               fn = 'annual_plant_n=18',
               num_species = 18, 
               dataset_complete = FALSE) # DEBUG
