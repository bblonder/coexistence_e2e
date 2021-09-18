library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(wesanderson)
library(tibble)
library(ggpubr)

dir.create('outputs_figures')

fn_outputs = dir('outputs_statistical',pattern="results*",full.names = TRUE)

df_all = rbindlist(lapply(1:length(fn_outputs), function(i)
{
  df_this = read.csv(fn_outputs[i])
  
  name_this = gsub("\\.csv","",gsub("results_","",basename(fn_outputs[i])))
  
  df_this$dataset_short = name_this
  
  df_this$dataset = sprintf("%s (%d total states)",name_this, 2^df_this$num_species[1] * df_this$num_replicates_in_data[1])
  
  return(df_this)
}))

df_all_stats = df_all %>% 
  select(dataset_short,num_species, num_replicates_in_data, num_cases) %>%
  #mutate(x.cutoff = 500 / (2^num_species * num_replicates_in_data)) %>%
  unique %>%
  arrange(num_species) %>%
  mutate(stochastic=dataset_short %in% c("sortie")) %>%
  mutate(empirical=dataset_short %in% c("cedar_creek"))

write.csv(df_all_stats,'outputs_figures/table_dataset_stats.csv',row.names=F)

 
make_plot <- function(yvar,xlab,ylab)
{
  ggplot(df_all,aes_string(x="frac",y=yvar,col="method")) +
    geom_point() +
    facet_wrap(~dataset) +
    theme_bw() +
    ylim(0,1) +
    scale_x_log10() +
    stat_summary(fun=mean, geom="line") +     
    ylab(xlab) +
    xlab(ylab) +
    theme(legend.position='bottom') +
    scale_color_manual(values=wes_palette("Darjeeling1")) +
    #geom_vline(mapping=aes(xintercept = x.cutoff), data=df_all_cutoff) +
    annotation_logticks(color='gray',alpha=0.5)
}

g_abundance = make_plot(yvar="abundance.r2", 
          xlab="Fraction of assemblages in training",
          ylab="R2 of abundance prediction")
ggsave(g_abundance, file='outputs_figures/g_abundance.png',width=7,height=5)

g_composition = make_plot(yvar="composition.balanced_accuracy", 
                        xlab="Fraction of assemblages in training",
                        ylab="Balanced accuracy of composition prediction")
ggsave(g_composition, file='outputs_figures/g_composition.png',width=7,height=5)

g_fs = make_plot(yvar="feasible.and.stable.balanced_accuracy", 
                          xlab="Fraction of assemblages in training",
                          ylab="Balanced accuracy of feasibility+stability prediction")
ggsave(g_fs, file='outputs_figures/g_fs.png',width=7,height=5)

g_richness = make_plot(yvar="richness.r2", 
                 xlab="Fraction of assemblages in training",
                 ylab="R2 of richness prediction")
ggsave(g_richness, file='outputs_figures/g_richness.png',width=7,height=5)






plot_obs_pred <- function(obs, pred, frac, limits)
{
  z = (as.matrix(obs) - as.matrix(pred)) / mean(as.numeric(as.matrix(obs)))
  row.names(z) = 1:nrow(z)
  print(min(z))
  print(max(z))
  #z = z[1:100,]
  
  z = z %>%
    as_tibble %>%
    rowid_to_column(var="X") %>%
    gather(key="Y", value="Z", -1) %>%
    mutate(Y = as.numeric(factor(Y)))
  
  g = ggplot(z, aes(X, Y, fill= Z)) + 
    theme_minimal() +
    theme(panel.border = element_rect(colour = "black", size = 1, fill=NA)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(legend.position='bottom') +
    geom_tile() +
    scale_fill_gradient2(low='red',mid='white',high='blue',midpoint = 0, name='(Obs-Pred)/Obs',limits=limits) +
    xlab('Assemblage') +
    ylab('Species') +
    ggtitle(bquote(beta ~ "=" ~ .(frac)))
  
  return(g)
}
  



plot_obs_pred_scatter <- function(obs, pred, frac)
{
  obs = obs %>% as.matrix %>% as.numeric
  pred = pred %>% as.matrix %>% as.numeric
  
  df = data.frame(obs=obs, pred=pred)
  
  minval = min(c(obs, pred))
  maxval = max(c(obs, pred))
  
  g = ggplot(data = df, aes(x = obs, y = pred)) +
    theme_bw() + 
    geom_hex(aes(fill = stat(log10(count)))) +
    scale_fill_gradient(low="lightgray",high="darkorchid") +
    xlab("Observed abundance") + 
    ylab("Predicted abundance") +
    stat_smooth(method='lm') +
    ggtitle(bquote(beta ~ "=" ~ .(frac))) +
    coord_fixed() +
    xlim(minval, maxval) + 
    ylim(minval, maxval)
  
  return(g)
}


data_abund_annual_0.001_obs = read.csv('outputs_statistical/table_fn=annual_plant_i=153_method=naive_rep=3_frac=0.001438_abundance_obs.csv')
data_abund_annual_0.001_pred = read.csv('outputs_statistical/table_fn=annual_plant_i=153_method=naive_rep=3_frac=0.001438_abundance_pred.csv')

data_abund_annual_0.01_obs = read.csv('outputs_statistical/table_fn=annual_plant_i=208_method=e2e_rep=3_frac=0.012743_abundance_obs.csv')
data_abund_annual_0.01_pred = read.csv('outputs_statistical/table_fn=annual_plant_i=208_method=e2e_rep=3_frac=0.012743_abundance_pred.csv')

data_abund_annual_0.05_obs = read.csv('outputs_statistical/table_fn=annual_plant_i=236_method=e2e_rep=1_frac=0.054556_abundance_obs.csv')
data_abund_annual_0.05_pred = read.csv('outputs_statistical/table_fn=annual_plant_i=236_method=e2e_rep=1_frac=0.054556_abundance_pred.csv')




limits_this = c(-50,50)
g_annualplant_heatmap = ggarrange(
  plot_obs_pred(data_abund_annual_0.001_obs, data_abund_annual_0.001_pred, 0.0014, limits=limits_this),
  plot_obs_pred(data_abund_annual_0.01_obs, data_abund_annual_0.01_pred, 0.013, limits=limits_this),
  plot_obs_pred(data_abund_annual_0.05_obs, data_abund_annual_0.05_pred, 0.054, limits=limits_this),
  nrow=1,ncol=3,
  common.legend=TRUE,
  legend='bottom')


g_annualplant_hexbin = ggarrange(
  plot_obs_pred_scatter(data_abund_annual_0.001_obs, data_abund_annual_0.001_pred, 0.0014),
  plot_obs_pred_scatter(data_abund_annual_0.01_obs, data_abund_annual_0.01_pred, 0.013),
  plot_obs_pred_scatter(data_abund_annual_0.05_obs, data_abund_annual_0.05_pred, 0.054),
  nrow=1,ncol=3,
  common.legend=TRUE,
  legend='bottom')


ggsave(ggarrange(g_annualplant_heatmap, g_annualplant_hexbin,
                 nrow=2,ncol=1,labels='auto'), 
       file='outputs_figures/g_annualplant_abundance.png',
       width=10,height=10,
       bg = 'white')


