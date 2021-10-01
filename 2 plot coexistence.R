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
  
  df_this$name = name_this
  
  #df_this$dataset = sprintf("%s (%d total states)",name_this, 2^df_this$num_species[1] * df_this$num_replicates_in_data[1])
  
  return(df_this)
}))

nice_names = c(fly_gut="Fly gut", 
               human_gut="Human gut",
               )

df_all_stats = df_all %>% 
  select(name,num_species, num_replicates_in_data, num_cases) %>%
  unique %>%
  arrange(num_species) %>%
  mutate(stochastic=name %in% c("sortie")) %>%
  mutate(empirical=name %in% c("cedar_creek","fly"))

write.csv(df_all_stats,'outputs_figures/table_dataset_stats.csv',row.names=F)

 
make_plot_performance <- function(data,yvar,ylab)
{
  ggplot(data,
         aes_string(x="frac",y=yvar,col="method")) +
    geom_point() +
    facet_grid(name~sampling_strategy,switch='x') +
    theme_bw() +
    ylim(0,1) +
    scale_x_log10() +
    stat_summary(fun=mean, geom="line") +     
    xlab(expression(paste(beta, " (training fraction out of ",q %*% 2^n,")"))) +
    ylab(ylab) +
    theme(legend.position='bottom') +
    scale_color_manual(values=wes_palette("Darjeeling1")) +
    #geom_vline(mapping=aes(xintercept = x.cutoff), data=df_all_cutoff) +
    annotation_logticks(color='gray',alpha=0.5)
}







g_abundance = make_plot_performance(data=df_all,
                        yvar="abundance.r2", 
                        ylab=expression(paste(R^2, " of abundance prediction")))
ggsave(g_abundance, file='outputs_figures/g_abundance.png',width=7,height=10)


g_composition = make_plot_performance(data=df_all,
                        yvar="composition.balanced_accuracy", 
                        ylab="Balanced accuracy of composition prediction")
ggsave(g_composition, file='outputs_figures/g_composition.png',width=7,height=5)

g_fs = make_plot_performance(data=df_all,
                      yvar="feasible.and.stable.balanced_accuracy",
                      ylab="Balanced accuracy of feasibility+stability prediction")
ggsave(g_fs, file='outputs_figures/g_fs.png',width=7,height=5)

g_richness = make_plot_performance(data=df_all,
                       yvar="richness.r2",
                       ylab="R2 of richness prediction")
ggsave(g_richness, file='outputs_figures/g_richness.png',width=7,height=5)










plot_obs_pred_scatter <- function(list_data)
{
  obs = read.csv(list_data$fn_obs) %>% 
    as.matrix %>% 
    as.numeric
  pred = read.csv(list_data$fn_pred) %>% as.matrix %>% as.numeric
  
  frac = list_data$best_frac
  
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


pick_datasets <- function(df, name, fraction, method, sampling_strategy)
{
  best_frac = df %>%
    filter(name==name & frac >= fraction & method==method & sampling_strategy==sampling_strategy) %>%
    arrange(frac) %>% slice_head()
  
  possible_files = dir('outputs_statistical',pattern='table*',full.names = TRUE)
  
  ids = grep(pattern=sprintf('fn=%s',name),possible_files)
  possible_files = possible_files[ids]
  
  ids = grep(pattern=sprintf('method=%s',method),possible_files)
  possible_files = possible_files[ids]
  
  ids = grep(pattern=sprintf('frac=%f',best_frac$frac[1]),possible_files)
  possible_files = possible_files[ids]
  
  ids = grep(pattern=sprintf('sampling_strategy=%s',sampling_strategy),possible_files)
  possible_files = possible_files[ids]
  
  print(possible_files)
  
  if(length(possible_files)==2)
  {
    file_obs = possible_files[grep("obs",possible_files)]
    file_pred = possible_files[grep("pred",possible_files)]
    
    return(list(fn_obs = file_obs,
             fn_pred = file_pred,
             best_frac = best_frac$frac[1]))
  }
  else
  {
    return(NULL)
  }
}


plot_scatter_dataset <- function(name)
{
  plots_all = lapply(c(0.001,0.01,0.05), function(i)
  { 
    l_this = pick_datasets(df_all, name=name, fraction=i, method='e2e', sampling_strategy='mixed')
    if (!is.null(l_this))
    {
      plot_obs_pred_scatter(l_this)
    }
    else
    {
      return(NULL)
    }
  })
  
  return(ggarrange(plotlist = plots_all,nrow=1,ncol=3,common.legend = TRUE,legend='right'))
}

perf_plots = lapply(df_all_stats$name, plot_scatter_dataset)

ggarrange(plotlist=perf_plots,nrow=4,ncol=1)





