library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(wesanderson)
library(tibble)
library(ggpubr)

try(dir.create('outputs_figures'))

fn_outputs = dir('outputs_statistical',pattern="results*",full.names = TRUE)

df_all = rbindlist(lapply(1:length(fn_outputs), function(i)
{
  df_this = read.csv(fn_outputs[i])
  
  name_this = gsub("\\.csv","",gsub("results_","",basename(fn_outputs[i])))
  
  df_this$name = name_this
  
  #df_this$dataset = sprintf("%s (%d total states)",name_this, 2^df_this$num_species[1] * df_this$num_replicates_in_data[1])
  
  return(df_this)
}))

nice_names = c(`annual_plant`="Annual plant",
  `cedar_creek_plants`="Cedar Creek",
  `fly_gut`="Fly gut", 
  `glv_simulated`="GLV simulated",
  `human_gut`="Human gut",
  `mouse_gut`="Mouse gut", 
  `sortie-nd_plants`="SORTIE-ND" 
               )

# add NA removal counts
source('quantile_trim.R')
fns = c(`cedar_creek_plants`='data_cedar_creek/cedar_creek_2018.csv', 
        `sortie-nd_plants`='data_sortie/data_sortie.csv',
        `human_gut`='data_glv/assemblages_H_12.csv',
        `mouse_gut`='data_glv/assemblages_M_11.csv',
        `glv_simulated`='data_glv/assemblages_glv_16.csv',
        `annual_plant`='data_annual_plant/assemblages_annual_plant_18.csv',
        `fly_gut`='data_fly/data_fly.csv')

row_counts = sapply(fns, function(x) {nrow(read.csv(x))})
row_counts_trimmed = sapply(fns, function(x) {
  q = quantile_trim(read.csv(x))
  counts = q %>% 
    select(contains("star")) %>%
    na.omit %>%
    nrow
  return(counts)
  })
num_nas = row_counts - row_counts_trimmed

# calculate stats
df_all_stats = df_all %>% 
  select(name,num_species, num_replicates_in_data) %>%
  mutate(nice_name = nice_names[name]) %>%
  unique %>%
  mutate(deterministic=name %in% c("mouse_gut","human_gut","glv_simulated")) %>%
  mutate(empirical=name %in% c("cedar_creek_plants","fly_gut")) %>%
  mutate(num_na = num_nas[name]) %>%
  mutate(num_states = row_counts[name]) %>%
  select(name, nice_name, deterministic, empirical, n=num_species, q=num_replicates_in_data, qx2n=num_states, num_na=num_na) %>%
  arrange(nice_name)


write.csv(df_all_stats %>% select(-name),'outputs_figures/table_dataset_stats.csv',row.names=F)

# get nice names
nn = df_all_stats$nice_name
names(nn) = df_all_stats$name









make_plot_performance <- function(data,yvar,ylab)
{
  data_ss = data %>% filter(num_train > 5) # cut off the single samples
  
  ggplot(data_ss,
         aes_string(x="frac",y=yvar,col="method")) +
    geom_point(alpha=0.5) +
    facet_grid(name~sampling_strategy,switch='x',
               labeller=labeller(name=nn)) +
    theme_bw() +
    ylim(0,1) +
    scale_x_log10(limits=c(1e-5,1e0),labels = function(x) format(x, scientific = TRUE)) +
    stat_summary(fun=mean, geom="line") +     
    xlab(expression(paste(beta, " (training fraction out of ",q %*% 2^n,")"))) +
    ylab(ylab) +
    theme(legend.position='bottom') +
    scale_color_manual(values=wes_palette("Darjeeling1")) +
    #geom_vline(mapping=aes(xintercept = x.cutoff), data=df_all_cutoff) +
    annotation_logticks(color='gray',alpha=0.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}







g_performance_abundance = make_plot_performance(data=df_all,
                        yvar="abundance.r2", 
                        ylab=expression(paste(R^2, " of abundance prediction")))
ggsave(g_performance_abundance, file='outputs_figures/g_performance_abundance.png',width=8,height=8)


g_performance_composition = make_plot_performance(data=df_all,
                        yvar="composition.balanced_accuracy", 
                        ylab="Balanced accuracy of composition prediction")
ggsave(g_performance_composition, file='outputs_figures/g_performance_composition.png',width=8,height=8)

g_performance_fs = make_plot_performance(data=df_all,
                      yvar="feasible.and.stable.balanced_accuracy",
                      ylab="Balanced accuracy of feasibility+stability prediction")
ggsave(g_performance_fs, file='outputs_figures/g_performance_fs.png',width=8,height=8)

g_performance_richness = make_plot_performance(data=df_all,
                       yvar="richness.r2",
                       ylab="R2 of richness prediction")
ggsave(g_performance_richness, file='outputs_figures/g_performance_richness.png',width=8,height=8)










plot_obs_pred_scatter <- function(list_data)
{
  obs = read.csv(list_data$fn_obs) %>% 
    as.matrix %>% 
    as.numeric
  pred = read.csv(list_data$fn_pred) %>% as.matrix %>% as.numeric
  
  frac = list_data$best_frac
  name = list_data$name
  
  df = data.frame(obs=obs, pred=pred)
  
  minval = min(c(obs, pred))
  maxval = max(c(obs, pred))
  
  g = ggplot(data = df, aes(x = obs, y = pred)) +
    theme_bw() + 
    geom_hex(aes(fill = stat(log10(count)))) +
    scale_fill_gradient(low="lightgray",high="darkorchid") +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    xlab("Observed abundance") + 
    ylab("Predicted abundance") +
    stat_smooth(method='lm') +
    ggtitle(nn[name]) +
    #ggtitle(bquote(beta ~ "=" ~ .(frac))) +
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
  
  ids = grep(pattern=sprintf('method=%s_',method),possible_files)
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
             best_frac = best_frac$frac[1],
             name = name))
  }
  else
  {
    return(NULL)
  }
}


plot_scatter_dataset <- function(name)
{
  l_this = pick_datasets(df_all, name=name, fraction=1e-2, method='e2e', sampling_strategy='mixed')
  if (!is.null(l_this))
  {
    p = plot_obs_pred_scatter(l_this)
  }
  else
  {
    return(NULL)
  }
  
  return(p)
  #return(ggarrange(plotlist = plots_all,nrow=1,ncol=length(plots_all),common.legend = TRUE,legend='right'))
}

perf_plots = lapply(df_all_stats$name, plot_scatter_dataset)

g_performance_scatter = ggarrange(plotlist=perf_plots,align='hv',common.legend = TRUE,legend='bottom')
ggsave(g_performance_scatter, file='outputs_figures/g_performance_scatter.png',width=10,height=10)




