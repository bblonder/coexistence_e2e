library(ggplot2)
library(reshape)
library(dplyr)
library(ggpubr)

source('quantile_trim.R')

plot_data <- function(data,name)
{
  # count # of species
  n_sp = data %>% select(contains("star")) %>% ncol
  
  # reorder data
  data = data %>% arrange(across(all_of(1:n_sp)))
  
  # remove outliers for all datasets
  data = quantile_trim(data)
  # plot abundance structure

  data_in = data %>%
    select(all_of(1:n_sp)) %>%
    mutate(row=1:nrow(.)) %>%
    melt(id.vars=c('row'))
  
  data_out = data %>%
    select(contains("star")) %>%
    mutate(row=1:nrow(.)) %>%
    melt(id.vars=c('row'))
  
  g_in = ggplot(data_in, aes(x=variable,y=row,fill=factor(value))) + 
    geom_raster() +
    scale_fill_manual(values=c('white','orange'),name='Experiment') +
    theme_bw() +
    xlab('Species') +
    ylab('Assemblage') +
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank()) +
    ggtitle(name)
  
  g_out = ggplot(data_out, aes(x=variable,y=row,fill=value)) + 
    geom_raster() +
    scale_fill_gradient2(name='Abundance') +
    theme_bw() +
    xlab('Species') +
    ylab('Assemblage') +
    theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  ggsave(ggarrange(g_in, g_out,nrow=1,ncol=2,align='hv'),file=sprintf('outputs_figures/g_experiment_%s.png',name),width=6,height=8)
}

try(dir.create('outputs_figures'))

data_assemblages_cedar_creek_18 = read.csv('data_cedar_creek/cedar_creek_2018.csv')
g_cedar_creek = plot_data(data_assemblages_cedar_creek_18,'Cedar Creek')

data_sortie_9_3 = read.csv('data_sortie/data_sortie.csv')
g_sortie = plot_data(data_sortie_9_3,'SORTIE')

data_assemblages_H_12 = read.csv('data_glv/assemblages_H_12.csv')
g_human_gut = plot_data(data_assemblages_H_12,'Human gut')

data_assemblages_M_11 = read.csv('data_glv/assemblages_M_11.csv')
g_mouse_gut = plot_data(data_assemblages_M_11,'Mouse gut')

data_assemblages_glv_16 = read.csv('data_glv/assemblages_glv_16.csv')
g_glv_simulated = plot_data(data_assemblages_glv_16,'GLV simulated')

data_annual_plant_18 = read.csv('data_annual_plant/assemblages_annual_plant_18.csv')
g_annual_plant = plot_data(data_annual_plant_18,'Annual plant')

data_fly_5 = read.csv('data_fly/data_fly.csv')
g_fly = plot_data(data_fly_5,'Fly gut')


