balanced_accuracy_casewise_mean <- function(pred, obs)
{
  stopifnot(nrow(pred) == nrow(obs))
  stopifnot(ncol(pred) == ncol(obs))
  
  balanced_accuracy_all = sapply(1:nrow(pred), function(i) {
    df_this_assemblage = data.frame(pred=factor(pred[i,,drop=TRUE],levels=c(FALSE,TRUE)), obs=factor(obs[i,,drop=TRUE],levels=c(FALSE,TRUE)))

    confusionMatrix(df_this_assemblage$pred, df_this_assemblage$obs)$byClass["Balanced Accuracy"]
  })
  
  balanced_accuracy_mean = mean(balanced_accuracy_all, na.rm=T)
  return(balanced_accuracy_mean)
}

mean_absolute_error_casewise_mean <- function(pred, obs)
{
  stopifnot(nrow(pred) == nrow(obs))
  stopifnot(ncol(pred) == ncol(obs))
  
  ae_all = sapply(1:nrow(pred), function(i) {
    df_this_assemblage = data.frame(pred=as.numeric(pred[i,,drop=TRUE]), obs=as.numeric(obs[i,,drop=TRUE]))
    
    # mean absolute error
    ae = mean_absolute_error(df_this_assemblage$pred, df_this_assemblage$obs)
    return(ae)
  })
  
  mae = mean(ae_all, na.rm=T)
  return(mae)
}

mean_absolute_error <- function(pred, obs)
{
  mae = mean(abs(as.numeric(pred) - as.numeric(obs)), na.rm=T)
  
  return(mae)
}

#data(iris)
#i1 = iris %>% filter(Species=="versicolor") %>% select(1:3) > 3
#i2 = iris %>% filter(Species=="setosa") %>% select(1:3) > 3

#ia = iris %>% filter(Species=="versicolor") %>% select(1:3)
#ib = iris %>% filter(Species=="setosa") %>% select(1:3)

#balanced_accuracy_casewise_mean(i1, i2)
#mean_absolute_error_casewise_mean(pred=ia, obs=ib)
