library(dplyr)
library(data.table)
library(reshape)
library(ggplot2)

files = dir('outputs',pattern = '*out',full.names = TRUE)

process_row <- function(i, plot=FALSE)
{
  df_this = read.table(files[i], sep = '\t', skip=5, header=TRUE)
  
  #df_this[nrow(df_this),grep("Sdl.Abs.Den.",names(df_this),fixed=TRUE)]
  n_final_this = tail(df_this[,grep("Sapl.Abs.Den.",names(df_this),fixed=TRUE)],1) + tail(df_this[,grep("Adult.Abs.Den.",names(df_this),fixed=TRUE)],1)
  
  species_names = sapply(strsplit(names(n_final_this),split="\\.\\."), tail, 1)
  
  n_last_50yr_this = tail(df_this[,grep("Sapl.Abs.Den.",names(df_this),fixed=TRUE)],10) + tail(df_this[,grep("Adult.Abs.Den.",names(df_this),fixed=TRUE)],10)
  
  n_cv_last_50yr_this = apply(n_last_50yr_this, 2, function(x) { cv = sd(x)/mean(x); cv[is.nan(cv)]=0; cv   })
  
  n_cv_last_50yr_mean_this = mean(n_cv_last_50yr_this)
  
  fn_this = gsub("_[0-9]\\.out","",gsub("GMF_","",basename(files[i])))
  composition_initial_this = as.numeric(strsplit(fn_this, '-')[[1]])
  
  df_action = data.frame(t(composition_initial_this))
  names(df_action) = paste(species_names,"action",sep=".")
  
  df_outcome = data.frame(t(as.numeric(n_final_this)))
  names(df_outcome) = paste(species_names,"outcome",sep=".")
  
  df_out = cbind(df_action, environment.initial=0, df_outcome)
  
  df_melted_adult = df_this %>% 
    select(Step, Adult.Abs.Den..ACRU:Adult.Abs.Den..QURU) %>% as.data.frame %>%
    melt(id=1)
  
  if (plot==TRUE)
  {
    g = ggplot(df_melted_adult,aes(x=Step,y=value,col=variable)) +
      geom_line() +
      ggtitle(files[1]) +
      theme_bw()
    
    ggsave(g, file=sprintf('~/Downloads/%s.pdf', basename(files[i])))
  }
  
  cat(sprintf('%f\n',i/length(files)))
  
  return(df_out)
}

rows_all = lapply(1:length(files), process_row, plot=FALSE)
df_final = rbindlist(rows_all)

write.csv(df_final, file='data_forest_trees.csv', row.names=FALSE)
