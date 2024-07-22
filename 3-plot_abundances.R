library(ggplot2)
library(reshape)
library(dplyr)
library(ggpubr)
library(data.table)
library(stringr)
library(tibble)

source('utils/quantile_trim.R')

plot_data <- function(data,name,trim=TRUE)
{
  # count # of species
  n_sp = data %>% select(contains("outcome")) %>% ncol
  
  # reorder data
  data = data %>% arrange(across(all_of(1:n_sp)))
  
  # flag quantile outliers
  if (trim==TRUE)
  {
    data = quantile_max_trim(data)
  }
  # plot abundance structure

  data_in = data %>%
    select(contains("action")) %>%
    mutate(row=1:nrow(.)) %>%
    as.data.frame %>% 
    reshape::melt(id.vars=c('row'))
  
  data_out = data %>%
    select(contains("outcome")) %>%
    mutate(row=1:nrow(.)) %>%
    as.data.frame %>% 
    reshape::melt(id.vars=c('row'))
  
  g_in = ggplot(data_in, aes(x=variable,y=row,fill=value)) + 
    geom_raster() +
    scale_fill_gradient(name='Experimental\naction',low='white',high='orange',na.value='red') +
    #scale_fill_manual(values=c('white','orange'),labels=c('Absent','Present'),name='Experimental\naction') +
    theme_bw() +
    xlab('Species') +
    ylab('Experiment') +
    theme(axis.text.x=element_blank()) +
    ggtitle(name) +
    scale_y_continuous(expand=c(0,0),breaks=range(data_in$row))
  
  g_out = ggplot(data_out, aes(x=variable,y=row,fill=value)) + 
    geom_raster() +
    scale_fill_gradient(name='Outcome\nabundance',low='white',high='blueviolet',na.value='red') +
    theme_bw() +
    xlab('Species') +
    ylab('') +
    theme(axis.text.x=element_blank()) +
    scale_y_continuous(expand=c(0,0),breaks=range(data_in$row))
  
  return(list(g_in, g_out))
}

try(dir.create('outputs/figures'))



fn_outputs = dir('outputs/statistical',pattern="result.*\\.csv",full.names = TRUE)

df_all_raw = rbindlist(lapply(1:length(fn_outputs), function(i)
{
  df_this = read.csv(fn_outputs[i])
  
  name_this = gsub("\\.csv","",gsub("results_","",basename(fn_outputs[i])))
  
  df_this$name = name_this
  
  return(df_this)
}))

truncate_name <- function(fn) {gsub('\\.csv','',gsub('outputs/statistical/results_','',fn))}

names_nice = truncate_name(fn_outputs)
names_nice_orig = names_nice
names_nice = str_to_sentence(gsub("_"," ",names_nice))
names_nice = gsub("Grassland ","Grassland\n",names_nice)
names(names_nice) = names_nice_orig

fns = sprintf('data/%s/data_%s.csv',names(names_nice),names(names_nice))
fns = gsub('data/human_gut','data/human_and_mouse_gut',fns)
fns = gsub('data/mouse_gut','data/human_and_mouse_gut',fns)
names(fns) = names_nice

lapply(1:length(fns), function(i) {
  print(i)
  data_this_raw = read.csv(fns[i])
    
  data_this = data_this_raw %>% 
    select(contains("initial"))
  
  # make name pairs
  names_to_insert = names(data_this)
  data_this = cbind(matrix(gsub("\\.initial","",names_to_insert),nrow=nrow(data_this),ncol=length(names_to_insert),byrow=TRUE),
                    data_this)
  data_this = data_this[,c(seq(1,ncol(data_this),by=2),seq(1,ncol(data_this),by=2)+1)]
  environment = apply(data_this, 1, paste, collapse=" ")
  # now reset without the extra columns
  data_this = data_this_raw %>% 
    mutate(environment = environment)
  
  print(table(data_this$environment))
  
  data_this_by_env = data_this %>% 
    group_by(environment) %>%
    group_split
  
  plots_this = lapply(1:length(data_this_by_env), function(j) {
    name_final_this = ifelse(length(unique(data_this$environment)) > 1,
                             paste(names(fns)[i],data_this_by_env[[j]]$environment[1],sep="\n"),
                             names(fns)[i])
    plots_this_env = plot_data(data_this_by_env[[j]] %>% select(-environment), name_final_this)
    return(ggarrange(plotlist=plots_this_env,nrow=1,ncol=2,align='hv'))
  })
  
  ggsave(ggarrange(plotlist=plots_this,ncol=1),
         file=sprintf('outputs/figures/g_experiment_%s.png',names(fns)[i]),width=6,height=4*length(data_this_by_env))
  })
