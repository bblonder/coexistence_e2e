library(ranger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(caret)
library(randomForestSRC)


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

predict_rf <- function(yvar=NULL, assemblages, rows_train, method, num_species=log2(nrow(assemblages)))
{
  rows_test = setdiff(1:nrow(assemblages), rows_train)
  
  data_for_rf_training = assemblages[rows_train,]
  
  if (method=='e2e + reshuffle')
  {
    # shuffle the x data
    data_for_rf_training[,names(data_for_rf_training)[1:num_species]] = data_for_rf_training[sample(1:nrow(data_for_rf_training)),names(data_for_rf_training)[1:num_species]]
  }
  
  if (is.null(yvar))
  {
    yvar = names(data_for_rf_training)[grep("star",names(data_for_rf_training))]
    
    # convert to presence/absence data only
    data_for_rf_training[,yvar] = data_for_rf_training[,yvar] > 0
    for (var in yvar)
    {
      data_for_rf_training[,var] = factor(data_for_rf_training[,var],levels=c(FALSE,TRUE))
    }
    
    if (method %in% c('e2e','e2e + reshuffle'))
    {
      formula_this_multivar = formula(sprintf("Multivar(%s)~%s",paste(yvar,collapse=", "),paste(names(data_for_rf_training)[1:num_species],collapse="+")))
      
      m_rf_multivar = rfsrc(formula=formula_this_multivar,
                            data = data_for_rf_training,
                            forest = TRUE)
      
      values_predicted_raw = predict(object=m_rf_multivar, 
                                 newdata=assemblages[rows_test,1:num_species])
      values_predicted = as.data.frame(sapply(values_predicted_raw$classOutput, function(x) {as.logical(x$class)}))
    }
    else if (method=='naive')
    {
      # use the input presence/absences
      values_predicted = assemblages[rows_test,1:num_species] > 0
      m_rf_multivar = NULL
    }
    
    values_observed = assemblages[rows_test,yvar] > 0
    
    #print(str(values_predicted))
    #print(str(values_observed))
    
    confusion = confusionMatrix(factor(as.numeric(as.matrix(values_predicted))),factor(as.numeric(as.matrix(values_observed))))
    
    balanced_accuracy = confusion$byClass["Balanced Accuracy"]

    return(list(model=m_rf_multivar,
                pred=values_predicted, 
                obs=values_observed,
                ba=balanced_accuracy))
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

do_prediction_graphs <- function(input_file,fn,num_species)
{
  data = input_file %>%
    mutate(feasible.and.stable = factor(stable & feasible))
  
  results_table = expand.grid(rep=1:3, method=c('e2e','naive','e2e + reshuffle'), frac=c(0.001,0.01,0.1,0.5,0.9))
  results_table$richness.r2 = NA
  results_table$feasible.and.stable.balanced_accuracy = NA
  results_table$composition.balanced_accuracy = NA
  
  for (i in 1:nrow(results_table))
  {
    message(i/nrow(results_table))
    rows_train = sample(x=1:nrow(data), size=ceiling(nrow(data)*results_table$frac[i]))
    
    prediction_richness = predict_rf(yvar = 'richness',
                                         assemblages = data,
                                         rows_train = rows_train,
                                         method = results_table$method[i],
                                         num_species=num_species)
    
    prediction_fs = predict_rf(yvar = 'feasible.and.stable',
                                   assemblages = data,
                                   rows_train = rows_train,
                                   method = results_table$method[i],
                                   num_species=num_species)
    
    prediction_composition = predict_rf(yvar = NULL,
                               assemblages = data,
                               rows_train = rows_train,
                               method = results_table$method[i],
                               num_species=num_species)
    
    results_table$richness.r2[i]=prediction_richness$r2
    results_table$feasible.and.stable.balanced_accuracy[i]=prediction_fs$ba
    results_table$composition.balanced_accuracy[i]=prediction_composition$ba
  }
  
  # MAKE PLOTS
  g_richness = ggplot(results_table, aes(x=factor(frac),y=richness.r2,col=method)) +
    geom_boxplot() +
    theme_bw() +
    ylim(0,1) +
    ylab("R2 for predicted richness") +
    xlab("Fraction of assemblages in training")
  ggsave(g_richness, file=sprintf('outputs/g_richness_%s.pdf',fn),width=8,height=6)
  
  g_fs = ggplot(results_table, aes(x=factor(frac),y=feasible.and.stable.balanced_accuracy,col=method)) +
    geom_boxplot() +
    theme_bw() +
    ylim(0,1) +
    ylab("Balanced accuracy of feasible/stable prediction") +
    xlab("Fraction of assemblages in training")
  ggsave(g_fs, file=sprintf('outputs/g_fs_%s.pdf',fn),width=8,height=6)
  
  
  g_pa = ggplot(results_table, aes(x=factor(frac),y=composition.balanced_accuracy,col=method)) +
    geom_boxplot() +
    theme_bw() +
    ylim(0,1) +
    ylab("Balanced accuracy of presence/absence prediction") +
    xlab("Fraction of assemblages in training")
  ggsave(g_pa, file=sprintf('outputs/g_pa_%s.pdf',fn),width=8,height=6)
}


# do analyses
data_assemblages_glv_10 = read.csv('data_glv/assemblages_glv_10.csv')
do_prediction_graphs(data_assemblages_glv_10,'glv_10',num_species = 10)


data_assemblages_H_12 = read.csv('data_glv/assemblages_H_12.csv')
do_prediction_graphs(data_assemblages_H_12,'H_12',num_species = 12)


data_assemblages_M_11 = read.csv('data_glv/assemblages_M_11.csv')
do_prediction_graphs(data_assemblages_M_11,'M_11',num_species = 11)

data_assemblages_CC_18 = read.csv('data_cedar_creek/cedar_creek_2018.csv')
do_prediction_graphs(data_assemblages_CC_18,'CC_18',num_species = 18)
154/2^18

