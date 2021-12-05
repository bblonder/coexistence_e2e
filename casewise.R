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

r2_casewise_mean <- function(pred, obs)
{
  stopifnot(nrow(pred) == nrow(obs))
  stopifnot(ncol(pred) == ncol(obs))
  
  r2_all = sapply(1:nrow(pred), function(i) {
    df_this_assemblage = data.frame(pred=as.numeric(pred[i,,drop=TRUE]), obs=as.numeric(obs[i,,drop=TRUE]))
    
    # avoid perfect fit warnings
    suppressWarnings(r2 <- summary(lm(pred~obs,data=df_this_assemblage))$r.squared)
  })
  
  r2_mean = mean(r2_all, na.rm=T)
  return(r2_mean)
}

#data(iris)
#i1 = iris %>% filter(Species=="versicolor") %>% select(1:3) > 3
#i2 = iris %>% filter(Species=="setosa") %>% select(1:3) > 3

#balanced_accuracy_casewise_mean(i1, i2)
