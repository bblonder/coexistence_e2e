library(dplyr)
library(data.table)

files = dir('outputs',pattern = '*out',full.names = TRUE)

process_row <- function(i)
{
  df_this = read.table(files[i], sep = '\t', skip=5, header=TRUE)
  
  #df_this[nrow(df_this),grep("Sdl.Abs.Den.",names(df_this),fixed=TRUE)]
  n_final_this = tail(df_this[,grep("Sapl.Abs.Den.",names(df_this),fixed=TRUE)],1) + tail(df_this[,grep("Adult.Abs.Den.",names(df_this),fixed=TRUE)],1)
  
  n_last_50yr_this = tail(df_this[,grep("Sapl.Abs.Den.",names(df_this),fixed=TRUE)],10) + tail(df_this[,grep("Adult.Abs.Den.",names(df_this),fixed=TRUE)],10)
  
  n_cv_last_50yr_this = apply(n_last_50yr_this, 2, function(x) { cv = sd(x)/mean(x); cv[is.nan(cv)]=0; cv   })
  
  n_cv_last_50yr_mean_this = mean(n_cv_last_50yr_this)
  
  fn_this = gsub("_[0-9]\\.out","",gsub("GMF_","",basename(files[i])))
  composition_initial_this = as.numeric(strsplit(fn_this, '-')[[1]])
  
  df1 = data.frame(t(composition_initial_this))
  names(df1) =letters[1:length(composition_initial_this)]
  
  df2 = data.frame(feasible=all(n_final_this>=0), 
                   stable=n_cv_last_50yr_mean_this<0.5,
                   richness = sum(n_final_this>=1))
  
  df3 = data.frame(t(as.numeric(n_final_this)))
  names(df3) = paste(letters[1:length(n_final_this)],"star",sep=".")
  
  df_out = cbind(df1, df2, df3)
  
  cat('.')
  
  return(df_out)
}

rows_all = lapply(1:length(files), process_row)
df_final = rbindlist(rows_all)

write.csv(df_final, file='data_sortie.csv', row.names=FALSE)
