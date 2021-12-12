library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(wesanderson)
library(tibble)
library(ggpubr)
library(MuMIn)
library(caret)


NUM_TRAIN_FOR_PLOTS = 62

if (!file.exists('outputs_figures'))
{
  dir.create('outputs_figures')
}

fn_outputs = dir('outputs_statistical',pattern="results*",full.names = TRUE)

df_all = rbindlist(lapply(1:length(fn_outputs), function(i)
{
  df_this = read.csv(fn_outputs[i])
  
  name_this = gsub("\\.csv","",gsub("results_","",basename(fn_outputs[i])))
  
  df_this$name = name_this

  return(df_this)
}))

nice_names = c(`annual_plant`="Annual plant",
  `cedar_creek_plants`="Cedar Creek",
  `fly_gut`="Fly gut", 
  `glv_simulated`="GLV simulated",
  `human_gut`="Human gut",
  `mouse_gut`="Mouse gut", 
  `sortie-nd_plants`="SORTIE-ND",
  `soil_bacteria`="Soil bacteria"
               )

# add NA removal counts
source('quantile_trim.R')
fns = c(`cedar_creek_plants`='data_cedar_creek/cedar_creek_2018.csv', 
        `sortie-nd_plants`='data_sortie/data_sortie.csv',
        `human_gut`='data_glv/assemblages_H_12.csv',
        `mouse_gut`='data_glv/assemblages_M_11.csv',
        `glv_simulated`='data_glv/assemblages_glv_16.csv',
        `annual_plant`='data_annual_plant/assemblages_annual_plant_18.csv',
        `fly_gut`='data_fly/data_fly.csv',
        `soil_bacteria`='data_friedman_gore/data_friedman_gore.csv')

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

# add some additional info
df_all = df_all %>% 
  mutate(nice_name = nice_names[name]) %>%
  mutate(empirical=name %in% c("cedar_creek_plants","fly_gut","soil_bacteria")) %>%
  mutate(num_na = num_nas[name]) %>%
  mutate(num_states = row_counts[name]) %>%
  mutate(richness_scaled_mae = richness_mae / num_species_dataset) %>%
  mutate(abundance_scaled_mae = abundance_mae_mean / abundance_q95_dataset)

# calculate stats
df_all_stats = df_all %>%
  select(nice_name, 
         n=num_species_dataset, 
         q=num_replicates_dataset, 
         num_states, 
         num_na, 
         empirical) %>%
  unique %>%
  arrange(nice_name)
write.csv(df_all_stats,'outputs_figures/table_dataset_stats.csv',row.names=F)











### PERFORMANCE OVERALL
varnames_nice = c(sampling_strategy='Sampling strategy',
                  num_species_dataset='Dataset, regional pool # species',
                  empirical='Dataset, from empirical study',
                  num_losses_mean_train='Experiment, # species lost',
                  abundance_final_skewness_mean_train='Outcome, abundance skewness',
                  nice_name='Dataset')

plot_visreg_perf <- function(yvar, ylab, ylim, num_train)
{
  # do regression analysis of performance
  df_perf_e2e = df_all %>% 
    filter(method=="e2e" & num_train==num_train) %>%
    mutate(empirical=factor(empirical)) %>%
    mutate(nice_name = factor(nice_name)) %>%
    mutate(sampling_strategy = factor(sampling_strategy)) %>%
    as.data.frame
  
  # remove NA responses
  df_perf_e2e = df_perf_e2e[!is.na(df_perf_e2e[,yvar]),]
  # clamp values above 1 to 1
  #df_perf_e2e[df_perf_e2e[,yvar] > 1,yvar] = 1
  
  xvars_all = c("empirical",
                "num_species_dataset", 
                "num_losses_mean_train",
                "abundance_final_skewness_mean_train", 
                "sampling_strategy",
                "nice_name")
  
  xvars_all_keep = xvars_all
  if(yvar=="feasible_and_stable_balanced_accuracy")
  {
    xvars_all_keep = setdiff(xvars_all,"empirical")
  }
  
  m_this = glm(formula = formula(sprintf("%s ~ %s", yvar, paste(xvars_all_keep,collapse=" + "))), 
                data=df_perf_e2e,family=quasibinomial)
  
  plots_all = vector(mode="list",length=length(xvars_all_keep))
  for (i in 1:length(xvars_all_keep))
  {
    plots_all[[i]] <- visreg(m_this, xvar=xvars_all_keep[i], gg=TRUE, overlay=FALSE, rug=FALSE, band=TRUE, scale='response') +
      xlab(varnames_nice[ xvars_all_keep[i] ]) +
      theme_bw() +
      theme(axis.text.x = element_text(
        angle = 30, hjust=1)) +
      ylab("") +
      ylim(ylim[1],ylim[2]) +
      labs(color='Dataset')
  }
  plots_arranged = ggarrange(plotlist=plots_all, common.legend = TRUE, legend='bottom',align='hv') +
    theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))
  plots_arranged = annotate_figure(plots_arranged, fig.lab=ylab)
  
  r2 = 1-(m_this$deviance/m_this$null.deviance)
  
  return(list(plots=plots_arranged, model=m_this, r2=r2, yvar=yvar, ncases=nrow(df_perf_e2e)))
}

vrp_richness = plot_visreg_perf(yvar="richness_scaled_mae",ylab="Scaled MAE of richness prediction (lower better)",ylim=c(0,0.5),num_train=NUM_TRAIN_FOR_PLOTS)
vrp_abundance = plot_visreg_perf("abundance_scaled_mae",ylab="Scaled MAE of abundance prediction (lower better)",ylim=c(0,0.5),num_train=NUM_TRAIN_FOR_PLOTS)
vrp_composition = plot_visreg_perf("composition_balanced_accuracy_mean",ylab="Balanced accuracy of composition prediction (higher better)",ylim=c(0,1),num_train=NUM_TRAIN_FOR_PLOTS)
vrp_feasiblestable = plot_visreg_perf("feasible_and_stable_balanced_accuracy",ylab="Balanced accuracy of feasibility+stability prediction (higher better)",ylim=c(0,1),num_train=NUM_TRAIN_FOR_PLOTS)

ggsave(vrp_richness$plots,file='outputs_figures/g_lm_richness.png',width=8,height=5)
ggsave(vrp_abundance$plots,file='outputs_figures/g_lm_abundance.png',width=8,height=5)
ggsave(vrp_composition$plots,file='outputs_figures/g_lm_composition.png',width=8,height=5)
ggsave(vrp_feasiblestable$plots,file='outputs_figures/g_lm_feasiblestable.png',width=8,height=5)

r2_lm = rbindlist(lapply(list(vrp_richness, vrp_abundance, vrp_composition, vrp_feasiblestable), function(x) { 
  data.frame(r2=x$r2, yvar=x$yvar) 
  }))

write.csv(r2_lm, file='outputs_figures/table_lm_r2.csv', row.names=FALSE)















make_plot_performance <- function(data,yvar,ylab,ylim)
{
  data_ss = data
  
  ggplot(data_ss,
         aes_string(x="num_train",y=yvar,col="method")) +
    geom_point(alpha=0.5) +
    facet_grid(name~sampling_strategy,switch='x',
               labeller=labeller(name=nice_names),scales='free_y') +
    theme_bw() +
    ylim(ylim[1],ylim[2]) +
    scale_x_log10(limits=c(1e1,1e4),labels = function(x) format(x, scientific = TRUE)) +
    stat_summary(fun=mean, geom="line") +     
    xlab("Number of training cases") +
    ylab(ylab) +
    theme(legend.position='bottom') +
    scale_color_manual(values=wes_palette("Darjeeling1")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}







g_performance_abundance = make_plot_performance(data=df_all,
                        yvar="abundance_scaled_mae", 
                        ylab="Scaled mean absolute error of abundance prediction (lower better)",
                        ylim=c(0,0.5))
ggsave(g_performance_abundance, file='outputs_figures/g_performance_abundance.png',width=8,height=9)


g_performance_composition = make_plot_performance(data=df_all,
                        yvar="composition_balanced_accuracy_mean", 
                        ylab="Mean balanced accuracy of composition prediction (higher better)",
                        ylim=c(0,1))
ggsave(g_performance_composition, file='outputs_figures/g_performance_composition.png',width=8,height=9)

g_performance_feasiblestable = make_plot_performance(data=df_all,
                      yvar="feasible_and_stable_balanced_accuracy",
                      ylab="Balanced accuracy of feasibility+stability prediction (higher better)",
                      ylim=c(0,1))
ggsave(g_performance_feasiblestable, file='outputs_figures/g_performance_feasiblestable.png',width=8,height=9)

g_performance_richness = make_plot_performance(data=df_all,
                       yvar="richness_scaled_mae",
                       ylab="Scaled mean absolute error of richness prediction (lower better)",
                       ylim=c(0,0.5))
ggsave(g_performance_richness, file='outputs_figures/g_performance_richness.png',width=8,height=9)






process_listdata <- function(listdata)
{
  obs = read.csv(listdata$fn_obs) %>% 
    as.matrix
  pred = read.csv(listdata$fn_pred) %>% 
    as.matrix

  stopifnot(nrow(obs)==nrow(pred))
  stopifnot(ncol(obs)==ncol(pred))
  
  rand_rows = sample(1:nrow(obs),min(100, nrow(obs)))
  obs = obs[rand_rows,,drop=FALSE]
  pred = pred[rand_rows,,drop=FALSE]
  
  matrix_all = matrix(data=NA,nrow=nrow(obs)*ncol(obs),ncol=3)
  matrix_all[,1] = as.numeric(obs)
  matrix_all[,2] = as.numeric(pred)
  matrix_all[,3] = sort(rep(1:nrow(obs),ncol(obs)))
  matrix_all = data.frame(matrix_all)
  names(matrix_all) = c('obs','pred','row')
  
  return(matrix_all)
}



plot_confusion_matrix <- function(listdata)
{
  matrix_all = process_listdata(listdata)
  
  xtab = table(matrix_all$pred, matrix_all$obs)
  
  cm = confusionMatrix(xtab)
  
  cm_melted = reshape::melt(as.matrix(cm),as.is=TRUE)
  names(cm_melted) = c("pred","obs","count")
  
  cm_melted$pred = factor(cm_melted$pred)
  cm_melted$obs = factor(cm_melted$obs)
  cm_melted$frac = cm_melted$count / sum(cm_melted$count)
  warning('make sure obs and pred are not swapped')
  
  g = ggplot(cm_melted,aes(x=pred,y=obs,fill=frac,label=round(frac,digits=2))) + 
    geom_tile() + geom_text() + 
    scale_fill_viridis_c(option='C',name='Fraction') +
    xlab("Predicted") +
    ylab("Observed") +
    theme_bw() +
    ggtitle(nice_names[listdata$name])
  
  return(g)
}

plot_obs_pred_scatter <- function(listdata, is_1d)
{
  matrix_all = process_listdata(listdata)
  
  allvals = c(as.numeric(matrix_all[,1]),as.numeric(matrix_all[,2]))
  
  minval = min(allvals,na.rm=TRUE)
  maxval = max(allvals,na.rm=TRUE)
  
  g = ggplot() +
    theme_bw() +
    scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    xlab("Observed") + 
    ylab("Predicted") +
    ggtitle(nice_names[listdata$name]) +
    coord_fixed() +
    xlim(minval, maxval) + 
    ylim(minval, maxval) +
    geom_point(data=matrix_all,aes(x=obs,y=pred),alpha=0.25,size=0.5) +
    scale_color_identity(guide = "legend")
  
  if (is_1d==FALSE)
  {
    g = g + geom_line(data=matrix_all,aes(x=obs,y=pred,group=row),stat='smooth',method='lm',alpha=0.5,color='olivedrab')
  }
  else
  {
    g = g + geom_line(data=matrix_all,aes(x=obs,y=pred),stat='smooth',method='lm',alpha=0.5,color='red')
  }  
  
  return(g)
}


pick_datasets <- function(df, name, num_train, method, sampling_strategy, response_variable)
{
  possible_files = dir('outputs_statistical',pattern='table*',full.names = TRUE)
  
  ids = grep(pattern=sprintf('fn=%s',name),possible_files)
  possible_files = possible_files[ids]
  
  ids = grep(pattern=sprintf('rep=%d_',1),possible_files)
  possible_files = possible_files[ids]
  
  ids = grep(pattern=sprintf('method=%s_',method),possible_files)
  possible_files = possible_files[ids]
  
  ids = grep(pattern=sprintf('num_train=%d',num_train),possible_files)
  possible_files = possible_files[ids]

  ids = grep(pattern=sprintf('sampling_strategy=%s',sampling_strategy),possible_files)
  possible_files = possible_files[ids]
  
  ids = grep(pattern=sprintf('%s',response_variable),possible_files)
  possible_files = possible_files[ids]
  
  if(length(possible_files)==2)
  {
    file_obs = possible_files[grep("obs",possible_files)]
    file_pred = possible_files[grep("pred",possible_files)]
    
    return(list(fn_obs = file_obs,
             fn_pred = file_pred,
             name = name))
  }
  else
  {
    return(NULL)
  }
}


plot_scatter_dataset <- function(name, response_variable, is_1d, is_continuous, num_train)
{
  l_this = pick_datasets(df_all, name=name, num_train=num_train, 
                               method='e2e', 
                               sampling_strategy='mixed', 
                               response_variable=response_variable)

  if (is_continuous==TRUE)
  {
    result = plot_obs_pred_scatter(l_this, is_1d=is_1d)
  } else
  {
    result = ggplot()
    try(result <- plot_confusion_matrix(l_this))
  }
  return(result)
}

perf_plots_abundance = lapply(unique(df_all$name), plot_scatter_dataset, response_variable='abundance', is_1d=FALSE, is_continuous=TRUE, num_train=NUM_TRAIN_FOR_PLOTS)
g_scatter_abundance = ggarrange(plotlist=perf_plots_abundance,align='hv',common.legend = TRUE,legend='bottom',nrow=2,ncol=4)
ggsave(g_scatter_abundance, file='outputs_figures/g_scatter_abundance.png',width=10,height=5)

perf_plots_richness = lapply(unique(df_all$name), plot_scatter_dataset, response_variable='richness', is_1d=TRUE, is_continuous=TRUE, num_train=NUM_TRAIN_FOR_PLOTS)
g_scatter_richness = ggarrange(plotlist=perf_plots_richness,align='hv',common.legend = TRUE,legend='bottom',nrow=2,ncol=4)
ggsave(g_scatter_richness, file='outputs_figures/g_scatter_richness.png',width=10,height=5)

perf_plots_composition = lapply(unique(df_all$name), plot_scatter_dataset, response_variable='composition', is_1d=TRUE, is_continuous=FALSE, num_train=NUM_TRAIN_FOR_PLOTS)
g_scatter_composition = ggarrange(plotlist=perf_plots_composition,align='hv',common.legend = TRUE,legend='bottom',nrow=2,ncol=4)
ggsave(g_scatter_composition, file='outputs_figures/g_scatter_composition.png',width=10,height=6)

perf_plots_feasiblestable = lapply(unique(df_all$name), plot_scatter_dataset, response_variable='fs', is_1d=TRUE, is_continuous=FALSE, num_train=NUM_TRAIN_FOR_PLOTS)
g_scatter_feasiblestable = ggarrange(plotlist=perf_plots_feasiblestable,align='hv',common.legend = TRUE,legend='bottom',nrow=2,ncol=4)
ggsave(g_scatter_feasiblestable, file='outputs_figures/g_scatter_feasiblestable.png',width=10,height=6)

