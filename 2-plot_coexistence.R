library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(wesanderson)
library(tibble)
library(ggpubr)
library(MuMIn)
library(caret)
library(visreg)
library(lme4)
library(vegan)
library(RColorBrewer)
library(MuMIn)
library(conflicted)
library(ggbiplot)
library(ggrepel)
library(car)
library(sjPlot)
library(ggheatmap)
library(DHARMa)
library(e1071)

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("mutate", "plyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("summarize", "dplyr")


colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###
if (!file.exists('outputs/figures'))
{
  dir.create('outputs/figures')  
}

###

fn_outputs = dir('outputs/statistical',pattern="result.*\\.csv",full.names = TRUE)

df_all_raw = rbindlist(lapply(1:length(fn_outputs), function(i)
{
  df_this = read.csv(fn_outputs[i])
  
  name_this = gsub("\\.csv","",gsub("results_","",basename(fn_outputs[i])))
  
  df_this$name = name_this
  
  return(df_this)
}))

names_nice = c(`annual_plant`="Annual plant",
               `cedar_creek_plants`="Cedar Creek",
               `fly_gut`="Fly gut", 
               #`glv_simulated`="GLV simulated",
               `human_gut`="Human gut",
               `mouse_gut`="Mouse gut", 
               `sortie-nd_plants`="SORTIE-ND",
               `soil_bacteria`="Soil bacteria"
)

varnames_nice = c(experimental_design='Experimental design',
                  num_species_dataset='Dataset, total # species',
                  type='Dataset, type',
                  num_losses_mean='Outcome, mean # species lost',
                  abundance_skewness_mean='Outcome, mean abundance skewness',
                  nice_name='Dataset'
                  )

experimental_design_nice = c(mixed='Mixed',
                             `low-2`='Doublets',
                             `low-3`='Doublets, triplets',
                             `high-1`='1-dropouts',
                             `high-2`='1-dropouts, 2-dropouts',
                             prior='Doublets, 1-dropouts; then mixed')

methods_nice = c(
                             rf='Random forest',
                             glv='GLV predictions',
                             glv_rf='Random forest on GLV residuals',
                             glv_rf_full='Random forest + GLV predictions',
                             sequential_rf='Random forest, sequential',
                             naive=' Naïve (mean abundance)')

# add NA removal counts
source('utils/quantile_trim.R')
fns = c(`cedar_creek_plants`='data/cedar_creek/cedar_creek_2018.csv', 
        `sortie-nd_plants`='data/sortie/data_sortie.csv',
        `human_gut`='data/glv/assemblages_H_12.csv',
        `mouse_gut`='data/glv/assemblages_M_11.csv',
        `annual_plant`='data/annual_plant/assemblages_annual_plant_18.csv',
        `fly_gut`='data/fly/data_fly.csv',
        `soil_bacteria`='data/friedman_gore/data_friedman_gore.csv')

row_counts = sapply(fns, function(x) {nrow(read.csv(x))})
row_counts_trimmed = sapply(fns, function(x) {
  data = read.csv(x)
  
  data = quantile_max_trim(data)
  counts = data %>% 
    select(contains("star")) %>%
    na.omit %>%
    nrow
  return(counts)
})
num_nas = row_counts - row_counts_trimmed

num_species = sapply(fns, function(x) {
  data = read.csv(x)
  
  n_sp = data %>% 
    select(contains("star")) %>% 
    ncol
  return(n_sp)
})

num_combos = sapply(fns, function(x) {
  data = read.csv(x)
  
  n_sp = data %>% 
    select(contains("star")) %>% 
    ncol
  
  num_combos = data %>%
    select(all_of(1:n_sp)) %>%
    unique %>% 
    nrow
  
  return(num_combos)
})

abundance_q95 = sapply(fns, function(x) {
  data = read.csv(x)
  
  q_95 = data %>% 
    select(contains("star")) %>% 
    as.matrix %>%
    as.numeric %>%
    quantile(0.95,na.rm=T) %>%
    as.numeric
  return(q_95)
})

# add some additional info
df_all = df_all_raw %>% 
  mutate(name_nice = names_nice[name]) %>%
  mutate(experimental_design_nice = factor(experimental_design_nice[experimental_design],levels=experimental_design_nice,ordered=TRUE)) %>%
  mutate(method_nice = factor(methods_nice[method],levels=methods_nice,ordered=TRUE)) %>%  
  
  mutate(empirical=name %in% c("cedar_creek_plants","fly_gut","soil_bacteria")) %>%
  mutate(num_na = num_nas[name]) %>%
  mutate(num_cases = row_counts[name]) %>%
  mutate(num_combos = num_combos[name]) %>%
  mutate(num_species_dataset = num_species[name]) %>%
  mutate(abundance_q95_dataset = abundance_q95[name]) %>%
  
  mutate(richness_mae_test_scaled = 
           richness_mae_test / num_species_dataset) %>%
  mutate(abundance_mae_mean_test_scaled = 
           abundance_mae_mean_test / abundance_q95_dataset) %>%
  mutate(one_minus_composition_balanced_accuracy_mean_test = 
           1 - composition_balanced_accuracy_mean_test) %>%
  mutate(abundance_mae_mean_test_scaled_clipped = 
           ifelse(abundance_mae_mean_test_scaled > 10, NA, abundance_mae_mean_test_scaled))


# calculate stats
df_all_stats = df_all %>%
  select(name_nice, 
         num_species_dataset,
         num_cases,
         num_combos,
         abundance_q95_dataset,
         num_na, 
         empirical) %>%
  mutate(possible_num_combos = 2^num_species_dataset) %>%
  mutate(abundance_q95_dataset = format(abundance_q95_dataset, nsmall=1, digits=1, scientific=FALSE)) %>%
  unique %>%
  arrange(tolower(name_nice))
write.csv(df_all_stats,'outputs/figures/table_dataset_stats.csv',row.names=F)






#### infrastructure to select cases
# because the outputs are too big, first download them
# cd /global/scratch/users/benjaminblonder
# find love/outputs/statistical | grep -e "method=rf\\_" | grep -e "num\\_train=89" | grep -e "experimental\\_design=mixed" | grep -e "response=abundance" | zip -@ test89.zip
# copy test.zip to local directory, move contents to outputs/statistical
# also find love/outputs/statistical | grep -e "method=rf\\_" | grep -e "num\\_train=264" | grep -e "experimental\\_design=mixed" | grep -e "response=abundance" | zip -@ test264.zip


source('utils/pick_datasets.R')
source('utils/log_seq.R')
possible_num_train = ceiling(log_seq(1e1,1e4,length.out=20))

cases_all = data.frame(name=names(names_nice), fn_assemblages=c('data/annual_plant/assemblages_annual_plant_18.csv',
                                                                'data/cedar_creek/cedar_creek_2018.csv',
                                                                'data/fly/data_fly.csv',
                                                                'data/glv/assemblages_H_12.csv',
                                                                'data/glv/assemblages_M_11.csv',
                                                                'data/sortie/data_sortie.csv',
                                                                'data/friedman_gore/data_friedman_gore.csv'))
cases = cases_all %>%
  filter(name %in% c('annual_plant','human_gut','mouse_gut','sortie-nd_plants')) # only do the complete datasets











g_abundance_mae_mean_test_scaled_clipped_by_num_train = ggplot(df_all %>%
         filter(experimental_design=='mixed'), 
       aes(x=num_train,
           y=abundance_mae_mean_test_scaled_clipped,
           color=name_nice,fill=name_nice)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  facet_wrap(~method_nice) +
  scale_color_brewer(palette='Set2',name='Dataset') +
  scale_fill_brewer(palette='Set2',name='Dataset') +
  geom_smooth(method='lm') +
  scale_x_log10() +
  scale_y_log10(breaks=c(0.01,0.05,0.2,1,10),limits=c(0.01,10)) +
  xlab("Number of training cases") +
  ylab("Scaled error") +
  theme(legend.position='right')

g_abundance_mae_mean_test_scaled_clipped_by_experimental_design = ggplot(df_all %>%
         filter(num_train==89 & method=='rf'), 
       aes(x=name_nice,
           y=abundance_mae_mean_test_scaled_clipped,
           color=experimental_design_nice,
           fill=experimental_design_nice,
           group=paste(experimental_design_nice,name_nice))) +
  #geom_boxplot(alpha=0.25,outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.5),alpha=0.5) +
  theme_bw() +
  scale_color_brewer(palette='Set1',name='Experimental design') +
  scale_fill_brewer(palette='Set1',name='Experimental design') +
  scale_y_log10(breaks=c(0.01,0.05,0.2,1,10),limits=c(0.01,10)) +
  xlab("Dataset") +
  ylab("Scaled error") +
  theme(legend.position='right')

g_abundance_mae_mean_test_scaled_clipped_by_method = ggplot(df_all %>%
                                                                           filter(num_train==89 & experimental_design=='mixed'), 
                                                                         aes(x=name_nice,
                                                                             y=abundance_mae_mean_test_scaled_clipped,
                                                                             color=method_nice,
                                                                             fill=method_nice,
                                                                             group=paste(method_nice,name_nice))) +
  #geom_boxplot(alpha=0.25,outlier.shape = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.5),alpha=0.75) +
  theme_bw() +
  scale_color_manual(values=brewer.pal(7,'Set3')[-2],name='Method') +
  scale_fill_manual(values=brewer.pal(7,'Set3')[-2],name='Method') +
  scale_y_log10(breaks=c(0.01,0.05,0.2,1,10),limits=c(0.01,10)) +
  xlab("Dataset") +
  ylab("Scaled error") +
  theme(legend.position='right')

g_abundance_mae = ggarrange(plotlist=list(
  g_abundance_mae_mean_test_scaled_clipped_by_num_train,
  g_abundance_mae_mean_test_scaled_clipped_by_method, 
  g_abundance_mae_mean_test_scaled_clipped_by_experimental_design),
  nrow=3,labels='auto')

ggsave(g_abundance_mae,
  width=8,height=9,file='outputs/figures/g_abundance_mae.pdf')
ggsave(g_abundance_mae,
       width=8,height=9,file='outputs/figures/g_abundance_mae.png')









### explanations
df_all_for_regression_dataset = df_all %>%
  filter(method=='rf' & num_train %in% c(21,89,264)) %>%
  mutate(type = factor(empirical,levels=c(FALSE,TRUE),labels = c('Simulated','Empirical')))




# 
# 
# g_dataset_properties_num_species = ggplot(df_all_for_regression_dataset, 
#        aes(color=factor(num_train), 
#            fill=factor(num_train),
#            y=abundance_mae_mean_test_scaled_clipped,
#            x=num_species_dataset)) +
#   #geom_boxplot(alpha=0.25,outlier.shape = NA) +
#   geom_smooth(method='lm') +
#   geom_point(position=position_jitterdodge(jitter.width = 0.5),alpha=0.5) +
#   theme_bw() +
#   scale_color_brewer(palette='Set2',name='Number of training cases') +
#   scale_fill_brewer(palette='Set2',name='Number of training cases') +
#   scale_y_log10(breaks=c(0.01,0.05,0.2,1,10),limits=c(0.01,10)) +
#   xlab(varnames_nice["num_species_dataset"]) +
#   ylab("Scaled error") +
#   theme(legend.position='right')
# 
# g_dataset_properties_dataset_type = ggplot(df_all_for_regression_dataset, 
#        aes(color=factor(num_train), 
#            fill=factor(num_train),
#            y=abundance_mae_mean_test_scaled_clipped,
#            x=type)) +
#   #geom_boxplot(alpha=0.25,outlier.shape = NA) +
#   geom_smooth(method='lm') +
#   geom_point(position=position_jitterdodge(jitter.width = 0.5),alpha=0.5) +
#   theme_bw() +
#   scale_color_brewer(palette='Set2',name='Number of training cases') +
#   scale_fill_brewer(palette='Set2',name='Number of training cases') +
#   scale_y_log10(breaks=c(0.01,0.05,0.2,1,10),limits=c(0.01,10)) +
#   xlab(varnames_nice["type"]) +
#   ylab("Scaled error") +
#   theme(legend.position='right')
# 
# g_dataset_properties_num_losses  = ggplot(df_all_for_regression_dataset, 
#        aes(color=factor(num_train), 
#            fill=factor(num_train),
#            y=abundance_mae_mean_test_scaled_clipped,
#            x=num_losses_mean)) +
#   #geom_boxplot(alpha=0.25,outlier.shape = NA) +
#   geom_smooth(method='lm') +
#   geom_point(position=position_jitterdodge(jitter.width = 0.5),alpha=0.5) +
#   theme_bw() +
#   scale_color_brewer(palette='Set2',name='Number of training cases') +
#   scale_fill_brewer(palette='Set2',name='Number of training cases') +
#   scale_y_log10(breaks=c(0.01,0.05,0.2,1,10),limits=c(0.01,10)) +
#   xlab(varnames_nice["num_losses_mean"]) +
#   ylab("Scaled error") +
#   theme(legend.position='right')
# 
# g_dataset_properties_abundance_final_skewness = ggplot(df_all_for_regression_dataset, 
#        aes(color=factor(num_train), 
#            fill=factor(num_train),
#            y=abundance_mae_mean_test_scaled_clipped,
#            x=abundance_final_skewness_mean)) +
#   #geom_boxplot(alpha=0.25,outlier.shape = NA) +
#   geom_smooth(method='lm') +
#   geom_point(position=position_jitterdodge(jitter.width = 0.5),alpha=0.5) +
#   theme_bw() +
#   scale_color_brewer(palette='Set2',name='Number of training cases') +
#   scale_fill_brewer(palette='Set2',name='Number of training cases') +
#   scale_y_log10(breaks=c(0.01,0.05,0.2,1,10),limits=c(0.01,10)) +
#   xlab(varnames_nice["abundance_final_skewness_mean"]) +
#   ylab("Scaled error") +
#   theme(legend.position='right')
# 
# 
# ggsave(ggarrange(plotlist=list(g_dataset_properties_num_species,
#                         g_dataset_properties_dataset_type,
#                         g_dataset_properties_num_losses,
#                         g_dataset_properties_abundance_final_skewness),
#           common.legend = TRUE,align='hv',labels = 'auto',legend='bottom'),
#        file='outputs/figures/g_dataset_properties.pdf',width=7,height=7)
# 



m_dataset_properties = lmer(abundance_mae_mean_test_scaled_clipped ~ (1|name) + 
                              log10(num_train)*num_species_dataset*type*num_losses_mean*abundance_skewness_mean,
     data=df_all_for_regression_dataset)
r.squaredGLMM(m_dataset_properties)
res_dataset_properties = simulateResiduals(m_dataset_properties,n=1000)

plot_model_properties <- function(xvar)
{
  visreg(m_dataset_properties,xvar=xvar,by='num_train',
         gg=TRUE,overlay=TRUE,jitter=TRUE,band=FALSE,
         type='conditional',line=list(alpha=0.5),point=list(alpha=0.5)) +
    theme_bw() +
    scale_color_brewer(palette='Set2',name='Number of training cases') +
    scale_fill_brewer(palette='Set2',name='Number of training cases') +
    ylim(0,0.22) +
    ylab('Scaled error') +
    xlab(varnames_nice[xvar])
}

g_dataset_properties_visreg = ggarrange(plotlist=list(plot_model_properties('num_species_dataset'),
                        plot_model_properties('num_losses_mean'),
                        plot_model_properties('abundance_skewness_mean'),
                        plot_model_properties('type')),
                        labels='auto',
          nrow=2,ncol=2,align='hv',legend='bottom',common.legend = TRUE)

ggsave(g_dataset_properties_visreg,
       file='outputs/figures/g_dataset_properties_visreg.pdf',width=7,height=7)
ggsave(g_dataset_properties_visreg,
       file='outputs/figures/g_dataset_properties_visreg.png',width=7,height=7)




























make_obs_pred_abundance_figure <- function(num_train, method_this='rf',experimental_design_this='mixed')
{
  predictions_abundance_all_for_scatter = rbindlist(lapply(1:nrow(cases_all), function(j) {
    cat('.')
    z = pick_datasets(name=cases_all$name[j], 
                      method=method_this,
                      num_train=num_train,
                      experimental_design = experimental_design_this,
                      response_var = 'abundance')
    
    #return(z)}))
    
    z_all = rbindlist(lapply(1:length(z$files_pred_test), function(i)
    {
      #print(paste(cases_all$name[j], i))
      if(length(z) > 0)
      {
        df_pred_test_this = read.csv(z$files_pred_test[i])
        df_obs_test_this = read.csv(z$files_obs_test[i])
        
        df_this = cbind(reshape2::melt(df_pred_test_this) %>%
                          rename(abundance_pred=value),
                        reshape2::melt(df_obs_test_this) %>%
                          rename(abundance_obs=value) %>%
                          select(-variable)) %>%
          mutate(rep=i) %>%
          mutate(variable=gsub("\\.star","",variable))
        
        df_this$name = cases_all$name[j]
        
        return(df_this)
      }
      else
      {
        return(NULL)
      }
    }))
    
    #
    
    return(z_all)
  }))
  
  # reduce # of cases to make the plotting work
  predictions_abundance_all_for_scatter_small = predictions_abundance_all_for_scatter %>%
    group_by(name, rep, variable) %>%
    sample_n(size=min(100,n())) %>%
    mutate(rep.factor=factor(rep))
  
  g_abundance_scatter = ggplot(predictions_abundance_all_for_scatter_small, aes(x=abundance_pred,
                                                                                y=abundance_obs,
                                                                                color=variable)) + 
    geom_line(stat="smooth",method = "lm", se=FALSE,alpha=0.5,linewidth=0.75, mapping=aes(group=paste(variable, rep))) +
    geom_point(alpha=0.5) +
    theme_bw() +
    geom_abline(slope=1,linewidth=2,alpha=0.5) +
    scale_color_hue(l=50,h=c(10,350)) +
    facet_wrap(~name,scales='free',labeller=as_labeller(names_nice),nrow=2) +
    labs(color='Species') +
    xlab('Abundance (predicted)') +
    ylab('Abundance (observed)')
  
  ggsave(g_abundance_scatter, file=sprintf('outputs/figures/g_abundance_scatter_%d.pdf', num_train),width=12,height=8)
  ggsave(g_abundance_scatter, file=sprintf('outputs/figures/g_abundance_scatter_%d.png', num_train),width=12,height=8)
}

make_obs_pred_abundance_figure(num_train=89, method='rf', experimental_design='mixed')
make_obs_pred_abundance_figure(num_train=264, method='rf', experimental_design='mixed')


















### PRIORITIZATION

get_datasets <- function(num_train)
{
  datasets_preds = lapply(1:nrow(cases), function(i) {
    cat(i)
    datasets_preds_this = pick_datasets(name=cases$name[i], 
                                        method='rf',
                                        num_train=num_train,
                                        experimental_design = 'mixed',
                                        response_var = 'abundance')
  })
  names(datasets_preds) = cases$name
  return(datasets_preds)
}

datasets_preds_89 = get_datasets(89)
datasets_preds_264 = get_datasets(264)




calculate_training_properties <- function(df_exp)
{
  outcomes_abundances = df_exp %>% select(contains('star'))
  num_species = ncol(outcomes_abundances)
  experiments = df_exp[,1:num_species]
  
  skewnesses = apply(outcomes_abundances, 1, skewness)
  
  richness_outcome = apply(outcomes_abundances, 1, function(x) { length(which(x>1e-6))   })
  richness_experiment = apply(experiments, 1, function(x) { length(which(x>1e-6))   })
  
  mean_skewness = mean(skewnesses)
  
  df_results = data.frame(skewnesses,richness_experiment,richness_outcome) %>%
    mutate(num_losses = richness_experiment - richness_outcome) %>%
    summarize(num_losses_mean = mean(num_losses), 
              abundance_skewness_mean = mean(skewnesses))
  
  return(df_results)
}

predictions_best <- function(type, datasets_preds, quantile_cutoff, fn_assemblages, options=NULL)
{
  if (is.null(datasets_preds)) { return(NULL) }
  cat('+')
  
  assemblages = read.csv(fn_assemblages)
  
  outcome_abundances_ALL = assemblages %>%
    select(contains("star"))
  outcome_abundances_ALL = quantile_max_trim(outcome_abundances_ALL)
  num_species = ncol(outcome_abundances_ALL)
  
  result = lapply (1:nrow(datasets_preds), function(i)
  {
    cat('.')
    
    # get the test set experiments
    df_experiments_train = read.csv(datasets_preds$files_exp_train[i])
    df_experiments_test = read.csv(datasets_preds$files_exp_test[i])
    experiments = rbind(df_experiments_train,df_experiments_test)[,1:num_species]
    
    # the predicted values for these experiments
    outcome_abundances_PREDICTED_train = read.csv(datasets_preds$files_pred_train[i])
    outcome_abundances_PREDICTED_test = read.csv(datasets_preds$files_pred_test[i])
    outcome_abundances_PREDICTED = rbind(outcome_abundances_PREDICTED_train, outcome_abundances_PREDICTED_test) 
    
    print(data.frame(nrow_all = nrow(outcome_abundances_ALL),
                     nrow_train = nrow(outcome_abundances_PREDICTED_train), 
                     nrow_test = nrow(outcome_abundances_PREDICTED_test),
                     nrow_predicted = nrow(outcome_abundances_PREDICTED)))
    
    if (type=='shannons_H')
    {
      shannons_H_PREDICTED = diversity(outcome_abundances_PREDICTED,index='shannon')
      #print(summary(outcome_abundances_ALL))
      shannons_H_ALL = diversity(outcome_abundances_ALL,index='shannon')
      
      #print(data.frame(length(shannons_H_PREDICTED), length(shannons_H_ALL)))
      #plot(shannons_H_PREDICTED, shannons_H_ALL)
      
      experiments_best_PREDICTED = experiments[shannons_H_PREDICTED >= quantile(shannons_H_ALL, quantile_cutoff, na.rm=T), 1:num_species] %>%
        na.omit
      experiments_best_ACTUAL = assemblages[shannons_H_ALL >= quantile(shannons_H_ALL, quantile_cutoff, na.rm=T), 1:num_species] %>%
        na.omit
      
      print(data.frame(num.best.predicted = nrow(experiments_best_PREDICTED), 
                       num.best.actual = nrow(experiments_best_ACTUAL)))
    }
    else if (type=='total_abundance')
    {
      total_abundance_PREDICTED = rowSums(outcome_abundances_PREDICTED)
      total_abundance_ALL = rowSums(outcome_abundances_ALL)
      
      # hist(total_abundance_ALL,breaks=100)
      # abline(v=quantile(total_abundance_ALL,quantile_cutoff,na.rm=T),col='red')
      # print(table(total_abundance_ALL > quantile(total_abundance_ALL, quantile_cutoff, na.rm=T)))
      # 
      experiments_best_PREDICTED = experiments[total_abundance_PREDICTED >= quantile(total_abundance_ALL, quantile_cutoff, na.rm=T), 1:num_species] %>%
        na.omit
      experiments_best_ACTUAL = assemblages[total_abundance_ALL >= quantile(total_abundance_ALL, quantile_cutoff, na.rm=T), 1:num_species] %>%
        na.omit  
    }
    else if (type=='remove_unwanted')
    {
      id_unwanted = options$id_unwanted
      outcome_abundance_unwanted_PREDICTED = outcome_abundances_PREDICTED[,id_unwanted]
      outcome_abundance_unwanted_ALL = outcome_abundances_ALL[,id_unwanted]
      #print(summary(outcome_abundance_unwanted_PREDICTED))
      #print(table((experiments[,id_unwanted]==1)))
      #print(summary(outcome_abundance_unwanted_ALL))
      #print(quantile(outcome_abundance_unwanted_ALL[outcome_abundance_unwanted_ALL > 0], quantile_cutoff, na.rm=T))
      
      experiments_best_PREDICTED = experiments[(experiments[,id_unwanted]==1) & 
                                                 (outcome_abundance_unwanted_PREDICTED <= quantile(outcome_abundance_unwanted_ALL, quantile_cutoff, na.rm=T)), 1:num_species] %>%
        na.omit
      
      
      experiments_best_ACTUAL = assemblages[(assemblages[,id_unwanted]==1) & 
                                              (outcome_abundance_unwanted_ALL <= quantile(outcome_abundance_unwanted_ALL, quantile_cutoff, na.rm=T)), 1:num_species] %>%
        na.omit  
      
      #print(c(nrow(experiments_best_PREDICTED),nrow(experiments_best_ACTUAL)))
    }
    else
    {
      stop('`type` not found')
    }
    
    ids_experiments_best_PREDICTED = apply(experiments_best_PREDICTED, 1, paste, collapse="*")
    ids_experiments_best_ACTUAL = apply(experiments_best_ACTUAL, 1, paste, collapse="*")
    ids_experiments_ALL = apply(assemblages[,1:num_species], 1, paste, collapse="*")
    
    flag_best_PREDICTED = ids_experiments_ALL %in% ids_experiments_best_PREDICTED
    flag_best_ACTUAL = ids_experiments_ALL %in% ids_experiments_best_ACTUAL
    
    confusion_matrix = confusionMatrix(data=factor(flag_best_PREDICTED,levels=c(FALSE,TRUE)), 
                                       reference=factor(flag_best_ACTUAL,levels=c(FALSE,TRUE)))
    
    confusion_matrix_stats = confusion_matrix$byClass %>% t %>% as.data.frame
    
    training_properties = calculate_training_properties(df_experiments_train)
    
    # also get abundances
    abundances_best_PREDICTED = outcome_abundances_PREDICTED[flag_best_PREDICTED,]
    abundances_best_ACTUAL = outcome_abundances_ALL[flag_best_ACTUAL,]
    
    return(list(name=datasets_preds$name[1],
                type=type,
                num_species=num_species,
                stats=confusion_matrix_stats,
                training_properties=training_properties,
                experiments_ALL=ids_experiments_ALL,
                #experiments_best_ACTUAL_matrix = experiments_best_ACTUAL,
                #experiments_best_PREDICTED_matrix = experiments_best_PREDICTED,
                #experiments_best_ACTUAL=ids_experiments_best_ACTUAL,
                #experiments_best_PREDICTED=ids_experiments_best_PREDICTED,
                abundances_best_PREDICTED=abundances_best_PREDICTED,
                #abundances_best_ACTUAL=abundances_best_ACTUAL,
                abundances_ALL = outcome_abundances_ALL,
                #abundances_best_ACTUAL_FOR_PREDICTIONS = outcome_abundances_ALL[flag_best_PREDICTED,],
                #abundances_best_ACTUAL_FOR_ACTUAL = outcome_abundances_ALL[flag_best_ACTUAL,],
                indices_best_ACTUAL = flag_best_ACTUAL,
                indices_best_PREDICTED = flag_best_PREDICTED
                  ))
  })
  cat('\n')
  return(result)
}


check_removal <- function(datasets_preds, dataset_name_this)
{
  num_species_this = read.csv(cases$fn_assemblages[cases$name==dataset_name_this]) %>% select(contains("star")) %>% ncol
  cat(dataset_name_this)
  stats_unwanted_this = lapply(1:num_species_this, function(id_unwanted_this) {
    cat('.')
    predictions_best_remove_unwanted = predictions_best(datasets_preds=datasets_preds[[dataset_name_this]],
                                                        type='remove_unwanted',
                                                        options=list(id_unwanted=id_unwanted_this),
                                                        quantile_cutoff=0.0,
                                                        fn_assemblages=cases$fn_assemblages[cases$name==dataset_name_this])
    
    result_this = rbindlist(lapply(predictions_best_remove_unwanted, function(x) { cbind(x$stats,x$training_properties) }))
    result_this$id_unwanted = id_unwanted_this
    result_this$name = dataset_name_this
    return(list(result_this=result_this,predictions_best_remove_unwanted=predictions_best_remove_unwanted))
  })
  
  result_df = rbindlist(lapply(stats_unwanted_this, function(x) { x$result_this }))
  result_best = lapply(stats_unwanted_this, function(x) { x$predictions_best_remove_unwanted })
  return(list(result_df=result_df, result_best=result_best))
}
predictions_best_removal_all_89 = lapply(cases$name, check_removal, datasets_preds=datasets_preds_89) # this runs slowly due to the annual plant dataset which is very large
predictions_best_removal_all_264 = lapply(cases$name, check_removal, datasets_preds=datasets_preds_264) # this runs slowly due to the annual plant dataset which is very large




predictions_best_shannons_H_89 = lapply(1:length(datasets_preds_89),
                                     FUN=function(i)
                                     {
                                       predictions_best(datasets_preds=datasets_preds_89[[i]],
                                                        type='shannons_H',
                                                        quantile_cutoff=0.95,
                                                        fn_assemblages=cases$fn_assemblages[i])
                                     })
predictions_best_shannons_H_264 = lapply(1:length(datasets_preds_264),
                                        FUN=function(i)
                                        {
                                          predictions_best(datasets_preds=datasets_preds_264[[i]],
                                                           type='shannons_H',
                                                           quantile_cutoff=0.95,
                                                           fn_assemblages=cases$fn_assemblages[i])
                                        })


predictions_best_total_abundance_89 = lapply(1:length(datasets_preds_89),
                                          FUN=function(i)
                                          {
                                            predictions_best(datasets_preds=datasets_preds_89[[i]],
                                                             type='total_abundance',
                                                             quantile_cutoff=0.95,
                                                             fn_assemblages=cases$fn_assemblages[i])
                                          })

predictions_best_total_abundance_264 = lapply(1:length(datasets_preds_264),
                                             FUN=function(i)
                                             {
                                               predictions_best(datasets_preds=datasets_preds_264[[i]],
                                                                type='total_abundance',
                                                                quantile_cutoff=0.95,
                                                                fn_assemblages=cases$fn_assemblages[i])
                                             })





# plot performance
make_removal_summary  <- function(results, num_train) {
  result = lapply(results, function(x) {x[[1]]}) %>% 
    rbindlist %>% 
    mutate(num_train=num_train) %>%
    mutate(name_nice=names_nice[name])
  return(result)
}
df_summary_removal_89 = make_removal_summary(predictions_best_removal_all_89, 89)
df_summary_removal_264 = make_removal_summary(predictions_best_removal_all_264, 264)

df_summary_removal_worst_89 = df_summary_removal_89 %>% 
  select(name, id_unwanted, num_train, Sensitivity, Specificity) %>% 
  melt(id.vars=c("name","id_unwanted","num_train")) %>%
  mutate(value=ifelse(is.na(value),0,value)) %>%
  mutate(id_unwanted = letters[id_unwanted]) %>%
  mutate(name_nice = names_nice[name]) %>%
  filter(variable=='Specificity') %>% 
  ungroup %>% 
  group_by(name_nice, id_unwanted, num_train, variable) %>% 
  summarize(mean_specificity=mean(value)) %>% 
  filter(mean_specificity < 0.2)
df_summary_removal_worst_264 = df_summary_removal_264 %>% 
  select(name, id_unwanted, num_train, Sensitivity, Specificity) %>% 
  melt(id.vars=c("name","id_unwanted","num_train")) %>%
  mutate(value=ifelse(is.na(value),0,value)) %>%
  mutate(id_unwanted = letters[id_unwanted]) %>%
  mutate(name_nice = names_nice[name]) %>%
  filter(variable=='Specificity') %>% 
  ungroup %>% 
  group_by(name_nice, id_unwanted, num_train, variable) %>% 
  summarize(mean_specificity=mean(value)) %>% 
  filter(mean_specificity < 0.2)

g_unwanted_all = ggplot(rbind(df_summary_removal_89,df_summary_removal_264) %>%
                          select(num_train, name_nice, id_unwanted,Sensitivity, Specificity) %>% 
                          melt(id.vars=c("num_train","name_nice","id_unwanted")), 
                        aes(x=factor(num_train),y=value,color=name_nice)) + 
  facet_wrap(~variable,scales='free_x',labeller=as_labeller(c(names_nice,Sensitivity='Sensitivity',Specificity='Specificity'))) +
  geom_point(position=position_jitterdodge(jitter.width = 0.5),alpha=0.5) +
  geom_label_repel(show.legend=FALSE, max.overlaps = 20, data=rbind(df_summary_unwanted_worst_89,df_summary_unwanted_worst_264), mapping=aes(x=factor(num_train),y=mean_specificity,color=name_nice,label=id_unwanted),position=position_jitterdodge(jitter.width = 0.5),alpha=0.5) +
  theme_bw() +
  ggtitle("Remove unwanted species") +
  xlab("Number of training cases") + ylab("Value (train + test)") +
  scale_color_manual(values=colorBlindBlack8,name='Dataset') +
  theme(legend.position='bottom')


get_classification_stats <- function(predictions_best, num_train)
{
  df_best_experiments_classification_statistics = lapply(predictions_best, 
                                                         FUN=function(predictions_this)
                                                         {
                                                           if (is.null(predictions_this)) { return(NULL) }
                                                           dataset_name = predictions_this[[1]]$name
                                                           
                                                           stats_this = rbindlist(lapply(predictions_this, function(x) { cbind(x$stats,x$training_properties) }))
                                                           stats_this$name = dataset_name
                                                           return(stats_this)
                                                         }) %>% 
    rbindlist %>%
    mutate(num_train=num_train)
  
  return(df_best_experiments_classification_statistics)
}

plot_best_experiments_classification_statistics <- function(predictions_best_list, num_train_list, title)
{
  stopifnot(length(predictions_best_list) == length(num_train_list))
  
  df_best_experiments_classification_statistics = rbindlist(lapply(1:length(num_train_list), function(i) {
    x = get_classification_stats(predictions_best_list[[i]], num_train_list[[i]]) %>% 
      select(name, Sensitivity, Specificity, num_train) %>%
      melt(id.vars=c('name','num_train'))
    
    #x$value[is.na(x$value)] = 0
    return(x)
    
    }))
            
  ggplot(df_best_experiments_classification_statistics %>% 
           mutate(name_nice=names_nice[name]), 
         aes(x=factor(num_train),y=value,color=name_nice)) +
    facet_wrap(~variable) +
    #geom_boxplot() +
    geom_point(position=position_jitterdodge(jitter.width = 0.25),alpha=0.5) +
    ylim(0,1) +
    theme_bw() +
    xlab("Number of training cases") +
    ylab("Value (train + test)") +
    scale_color_manual(values=colorBlindBlack8,name='Dataset') +
    ggtitle(title)
}


g_best_experiments_classification_statistics_shannons_H = plot_best_experiments_classification_statistics(list(predictions_best_shannons_H_89,predictions_best_shannons_H_264),c(89,264), 
                                                                                                     ">95% quantile Shannon's H")
g_best_experiments_classification_statistics_total_abundance = plot_best_experiments_classification_statistics(list(predictions_best_total_abundance_89, predictions_best_total_abundance_264),c(89,264), 
                                                                                                          ">95% quantile total abundance")


g_best_experiments_classification_statistics = ggarrange(g_unwanted_all,
                                                          g_best_experiments_classification_statistics_shannons_H, 
                                                          g_best_experiments_classification_statistics_total_abundance,
                    nrow=2,ncol=2,labels='auto',common.legend = TRUE, legend='bottom')

ggsave(g_best_experiments_classification_statistics,
       file='outputs/figures/g_best_experiments_classification_statistics.pdf',
       width=9,height=9)
ggsave(g_best_experiments_classification_statistics,
       file='outputs/figures/g_best_experiments_classification_statistics.png',
       width=9,height=9)





# do a post-hoc model
df_prioritization_stats_removal = rbind(df_summary_removal_89,
    df_summary_removal_264) %>%
  mutate(problem='removal') %>%
  #mutate(num_train = factor(num_train)) %>%
  mutate(num_species_dataset = num_species[name]) %>%
  mutate(name_nice=names_nice[name]) %>%
  mutate(id_wanted=paste(name,id_unwanted)) # make unique


df_prioritization_stats_shannons_h = rbind(get_classification_stats(predictions_best_shannons_H_89, 89),
    get_classification_stats(predictions_best_shannons_H_264, 264)) %>%
  mutate(problem='shannons_h') %>%
  #mutate(num_train = factor(num_train)) %>%
  mutate(num_species_dataset = num_species[name]) %>%
  mutate(name_nice=names_nice[name])

df_prioritization_stats_abundance = rbind(get_classification_stats(predictions_best_total_abundance_89, 89),
                                           get_classification_stats(predictions_best_total_abundance_264, 264)) %>%
  mutate(problem='abundance') %>%
  #mutate(num_train = factor(num_train)) %>%
  mutate(num_species_dataset = num_species[name]) %>%
  mutate(name_nice=names_nice[name])




m_specificity_removal = lmer(Specificity ~ log10(num_train) + num_losses_mean + abundance_skewness_mean + num_species_dataset + (1|name_nice) + (1|id_unwanted),
     data=df_prioritization_stats_removal)
Anova(m_specificity_removal)
r.squaredGLMM(m_specificity_removal)
res_specificity_removal = simulateResiduals(m_specificity_removal,n=1000)
plot(res_specificity_removal)


m_specificity_shannons_h = lmer(Specificity ~ log10(num_train) + num_losses_mean + abundance_skewness_mean + num_species_dataset + (1|name_nice),
                                data=df_prioritization_stats_shannons_h)
Anova(m_specificity_shannons_h)
r.squaredGLMM(m_specificity_shannons_h)
res_specificity_shannons_h = simulateResiduals(m_specificity_shannons_h,n=1000)
plot(res_specificity_shannons_h)

m_specificity_abundance = lmer(Specificity ~ log10(num_train) + num_losses_mean + abundance_skewness_mean + num_species_dataset + (1|name_nice),
                                data=df_prioritization_stats_abundance)
Anova(m_specificity_abundance)
r.squaredGLMM(m_specificity_abundance)
res_specificity_abundance = simulateResiduals(m_specificity_abundance,n=1000)
plot(res_specificity_abundance)

plot_specificity_model_properties <- function(model_this, xvar)
{
  visreg(model_this,xvar=xvar,by='num_train',
         gg=TRUE,overlay=TRUE,jitter=TRUE,band=FALSE,
         type='conditional',line=list(alpha=0.5),point=list(alpha=0.5)) +
    theme_bw() +
    scale_color_brewer(palette='Set1',name='Number of training cases') +
    scale_fill_brewer(palette='Set1',name='Number of training cases') +
    ylab('Specificity') +
    xlab(varnames_nice[xvar]) +
    scale_y_continuous(breaks=c(0,0.5,1),limits=c(-0.2,1.2))
}

pl1 = lapply(c('num_species_dataset','num_losses_mean','abundance_skewness_mean'), 
             plot_specificity_model_properties, 
             model_this=m_specificity_removal)
pl1[[1]] = pl1[[1]] + ggtitle('Remove unwanted species')

pl2 = lapply(c('num_species_dataset','num_losses_mean','abundance_skewness_mean'), 
             plot_specificity_model_properties, 
             model_this=m_specificity_shannons_h)
pl2[[1]] = pl2[[1]] + ggtitle(">95% quantile Shannon's H")

pl3 = lapply(c('num_species_dataset','num_losses_mean','abundance_skewness_mean'), 
             plot_specificity_model_properties, 
             model_this=m_specificity_removal)
pl3[[1]] = pl3[[1]] + ggtitle(">95% quantile abundance")

g_specificity_models_visreg = ggarrange(plotlist=c(pl1, pl2, pl3),
          nrow=3,ncol=3,labels='auto',common.legend = TRUE, legend='bottom',align='hv')
ggsave(g_specificity_models_visreg, file='outputs/figures/g_specificity_models_visreg.pdf',width=9.5,height=9.5)
ggsave(g_specificity_models_visreg, file='outputs/figures/g_specificity_models_visreg.png',width=9.5,height=9.5)

# make plots of best experiments (old)
# 
# plot_best_experiments_error_rates <- function(df, extra_text="")
# {
#   df_short = df %>% 
#     group_split(rep) %>% 
#     #lapply(function(x) {x$actual}) %>%
#     lapply(function(x) {x$classification}) %>% 
#     as.data.frame
#   names(df_short) = 1:ncol(df_short)
#   df_long = reshape2::melt(as.matrix(df_short)) %>%
#     mutate(Var2=as.character(Var2)) %>%
#     #mutate(value=factor(value))
#     mutate(value=factor(value,levels=levels(df$classification))) %>%
#     mutate(Var2=factor(Var2,levels=1:10,labels=paste(1:10),ordered=TRUE))
#   
#   g = ggplot(df_long,aes(x=Var1,y=Var2,fill=value)) + 
#     geom_tile() + 
#     scale_fill_manual(values=c('red','orange3','black','blue'),name='Classification',drop=FALSE) + 
#     theme_bw() +
#     xlab('Experiment') +
#     ylab('Replicate') +
#     ggtitle(paste(df$name[1], df$type[1])) +
#     scale_x_continuous(expand=c(0,0),breaks=range(df_long$Var1)) +
#     ggtitle(paste(names_nice[df$name[1]], extra_text)) + 
#     theme(plot.title = element_text(size = 8)) +
#     theme(axis.text.x=element_text(size=3))
#   
#   return(g)
# }
# 
# plots_best_experiments_error_rates_shannons_H = lapply(1:length(best_experiments_matrix_shannons_H), function(i) {
#   plot_best_experiments_error_rates(best_experiments_matrix_shannons_H[[i]], 
#                                     extra_text = " - find >95% quantile Shannon's H")
# })
# 
# plots_best_experiments_error_rates_total_abundance = lapply(1:length(best_experiments_matrix_total_abundance), function(i) {
#   plot_best_experiments_error_rates(best_experiments_matrix_total_abundance[[i]], 
#                                     extra_text = " - find >95% quantile total abundance")
# })
# 
# plots_best_experiments_error_rates_remove_unwanted = lapply(1:length(best_experiments_matrix_remove_unwanted), function(i) {
#   plot_best_experiments_error_rates(best_experiments_matrix_remove_unwanted[[i]], 
#                                     extra_text = sprintf(" - remove %s", letters[ best_rows_unwanted$id_unwanted[i] ]))
# })
# 
# g_best_experiments_error_rates_all = ggarrange(plotlist=c(plots_best_experiments_error_rates_remove_unwanted,
#                                                           plots_best_experiments_error_rates_shannons_H,
#                                                           plots_best_experiments_error_rates_total_abundance
#   ),
#     nrow=3,ncol=4,
#     common.legend = TRUE,
#     labels=c('(a)',rep('',3),'(b)',rep('',3),'(c)',rep('',3)),
#     legend='bottom',
#     align='hv')
# ggsave(g_best_experiments_error_rates_all,
#        file='outputs/figures/g_best_experiments_error_rates_all.png',
#        width=12,height=8,dpi=1600)
# 
# 
# 
# 
# 
# plot_best_experiments_presence_absence <- function(fn_assemblages, df, extra_text="")
# {
#   cat("*")
#   
#   experiments_all = read.csv(fn_assemblages) 
#   num_species = ncol(experiments_all %>% select(contains("star")))
#   experiments_all = experiments_all[,1:num_species]
#   ids_experiments_all = apply(experiments_all, 1, paste, collapse="*")
#   
#   experiments_to_plot = experiments_all
#   # convert to character format
#   for (j in 1:ncol(experiments_to_plot))
#   {
#     experiments_to_plot[,j] = as.character(experiments_to_plot[,j])
#   }
#   
#   # pick the best output
#   df_count_pos = df %>% 
#     group_by(rep) %>% 
#     summarize(tp=sum(classification=='true positive')) %>%
#     ungroup
#   
#   id_best = df_count_pos$rep[which.max(df_count_pos$tp)]
#   #print(df_count_pos)
#   
#   df_this = df %>% filter(rep==id_best)
#   
#   # make sure the ordering of rows matches
#   stopifnot(all((df %>% filter(rep==id_best) %>% pull(id)) == ids_experiments_all))
#   
#   experiments_to_plot_final = matrix(rep(as.character(df_this$classification),ncol(experiments_to_plot)),
#                                      nrow=nrow(experiments_to_plot),ncol=ncol(experiments_to_plot))
#   experiments_to_plot_final[experiments_to_plot=="0"] = NA
#   
#   experiments_to_plot_long = reshape2::melt(as.matrix(experiments_to_plot_final)) %>%
#     mutate(value=factor(value,levels=levels(df$classification))) %>%
#     mutate(Var2=factor(Var2,levels=1:num_species,labels=letters[1:num_species]))
#   
#   g = ggplot(experiments_to_plot_long,aes(x=Var1,y=Var2,fill=value)) + 
#     geom_tile() + 
#     scale_fill_manual(values=c('red','orange3','black','blue'),name='Classification',drop=FALSE) +
#     theme_bw() +
#     xlab('Experiment') +
#     ylab('Species') +
#     scale_x_continuous(expand=c(0,0),breaks=range(experiments_to_plot_long$Var1)) +
#     ggtitle(paste(names_nice[df$name[1]], "-", "replicate", id_best, "-", extra_text)) + 
#     theme(plot.title = element_text(size = 8)) +
#     theme(axis.text.x=element_text(size=3))
#   
#   cat('\n')
#   
#   return(g)
# }
# 
# plots_best_experiments_presence_absence_shannons_H = lapply(1:length(best_experiments_matrix_shannons_H),
#                                                             FUN=function(i)
#                                                             {
#                                                               plot_best_experiments_presence_absence(
#                                                                 fn_assemblages=cases$fn_assemblages[i],
#                                                                 df=best_experiments_matrix_shannons_H[[i]],
#                                                                 extra_text = "\nfind >95% quantile Shannon's H")
#                                                             })
# 
# plots_best_experiments_presence_absence_total_abundance = lapply(1:length(best_experiments_matrix_total_abundance),
#                                                                  FUN=function(i)
#                                                                  {
#                                                                    plot_best_experiments_presence_absence(
#                                                                      fn_assemblages=cases$fn_assemblages[i],
#                                                                      df=best_experiments_matrix_total_abundance[[i]],
#                                                                      extra_text = "\nfind >95% quantile total abundance")
#                                                                  })
# 
# plots_best_experiments_presence_absence_remove_unwanted = lapply(1:length(best_experiments_matrix_remove_unwanted),
#                                                                  FUN=function(i)
#                                                                  {
#                                                                    plot_best_experiments_presence_absence(
#                                                                      fn_assemblages=cases$fn_assemblages[i],
#                                                                      df=best_experiments_matrix_remove_unwanted[[i]],
#                                                                      extra_text = sprintf("\nremove %s", letters[ best_rows_unwanted$id_unwanted[i] ]))
#                                                                  })
# 
# 
# g_best_experiments_presence_absence_all = ggarrange(plotlist=c(plots_best_experiments_presence_absence_remove_unwanted,
#                                                                plots_best_experiments_presence_absence_shannons_H,
#                                                                plots_best_experiments_presence_absence_total_abundance
#   ),
#     nrow=3,ncol=4,
#     common.legend = TRUE,
#     
#     labels=c('(a)',rep('',3),'(b)',rep('',3),'(c)',rep('',3)),
#     legend='bottom',
#     align='hv'
#     )
# 
# ggsave(g_best_experiments_presence_absence_all,
#        file='outputs/figures/g_best_experiments_presence_absence_all.png',
#        width=12,height=8,dpi=1600)
# 
# 









### SHOW matrices of best experiments

plot_best_experiments_new <- function(predictions_list, num_experiments_max=500)
{
  cat('.')
  predictions_this = predictions_list$predictions
  
  experiments_matrix = rbindlist(lapply(predictions_this$experiments_ALL, function(x) { 
    x = strsplit(x, split="\\*"); as.data.frame(t(as.numeric(unlist(x)))) }) )
  
  # do some resampling due to memory issues
  if (length(which(predictions_this$indices_best_ACTUAL==TRUE)) > num_experiments_max)
  {
    predictions_this$indices_best_ACTUAL = sample(which(predictions_this$indices_best_ACTUAL==TRUE), num_experiments_max, replace=FALSE)
  }
  if (length(which(predictions_this$indices_best_PREDICTED==TRUE)) > num_experiments_max)
  {
    predictions_this$indices_best_PREDICTED = sample(which(predictions_this$indices_best_PREDICTED==TRUE), num_experiments_max, replace=FALSE)
  }
  
  m_experiments_best_actual = experiments_matrix[predictions_this$indices_best_ACTUAL,]
  m_experiments_best_predicted = experiments_matrix[predictions_this$indices_best_PREDICTED,]
  
  m_outcomes_best_actual = predictions_this$abundances_ALL[predictions_this$indices_best_ACTUAL,]
  m_outcomes_best_predicted = predictions_this$abundances_ALL[predictions_this$indices_best_PREDICTED,]
  #m_outcomes_best_predicted = predictions_this$abundances_best_PREDICTED[indices_predicted_to_sample,] # assumes ordering was set up correctly in above code... which it is! but not an elegant solution
  
  cat('.')
  if (!is.null(m_outcomes_best_actual))
  {
    dist_actual = dist(m_outcomes_best_actual)
    order_actual = hclust(dist_actual)$order
  }
  else
  {
    order_actual = NULL
  }
  
  if (!is.null(m_outcomes_best_predicted) & !is.null(predictions_this))
  {
    if (predictions_this$stats$Specificity>0)
    {
      dist_predicted = dist(m_outcomes_best_predicted)
      order_predicted = hclust(dist_predicted)$order
    }
    else
    {
      order_predicted = NULL
    }
  }
  else
  {
    order_predicted = NULL
  }
  
  plot_matrix_heatmap <- function(m, ordering, type)
  {
    if (nrow(m) == 0) 
    { 
      return(ggplot()) 
    }
    else
    {
      names(m) = letters[1:ncol(m)]
      row.names(m) <- 1:nrow(m)
      
      if(!is.null(ordering))
      {
        m = m[ordering,]
      }
      
      m_melted = m %>% # Data wrangling
        rowid_to_column(var="Y") %>%
        gather(key="X", value="Z", -1)# %>%
      # Change Y to numeric
      #mutate(Y=as.numeric(gsub("V","",Y)))
      
      if (type=='experiment')
      {
        m_melted$Z = factor(m_melted$Z,levels=c(0,1))
      }
      
      g = ggplot(m_melted, aes(x=X,y=Y,fill=Z)) +
        geom_tile() +
        theme_bw() + 
        scale_x_discrete(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) +
        theme(axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.title.y=element_blank(),
              axis.title.x=element_blank(),
              plot.margin=margin(t=20,l=10))
      
      if (type=='experiment')
      {
        g = g +
          scale_fill_manual(values=c('lightgray','darkorange4'),name='Experimental action')
      }
      else if (type=='outcome')
      {
        g = g +
          scale_fill_gradientn(colors=brewer.pal(9,"YlGnBu"),
                               name='Outcome abundance',trans='log',
                               labels=function(x) sprintf("%.1f", x),
                               na.value = 'lightgray')
      }
      
      return(g)
    }
  }
  
  cat('.')
  if (!is.null(m_experiments_best_actual))
  {
    g_e_a = plot_matrix_heatmap(m_experiments_best_actual,order_actual,type='experiment') + ggtitle('Actual best experiments')
  }
  else
  {
    g_e_a = ggplot() + ggtitle('Actual best experiments')
  }
  if (!is.null(m_experiments_best_predicted))
  {
    g_e_p = plot_matrix_heatmap(m_experiments_best_predicted,order_predicted,type='experiment') + ggtitle('Predicted best experiments')
  }
  else
  {
    g_e_p = ggplot() + ggtitle('Predicted best experiments')
  }
  print(str(m_outcomes_best_actual))
  if (!is.null(m_outcomes_best_actual))
  {
    g_o_a = plot_matrix_heatmap(m_outcomes_best_actual,order_actual,type='outcome') + ggtitle('Actual outcome of actual best experiments')
  }
  else
  {
    g_o_a = ggplot() + ggtitle('Actual outcome of actual best experiments')
  }
  if (!is.null(m_outcomes_best_predicted))
  {
    g_o_p = plot_matrix_heatmap(m_outcomes_best_predicted,order_predicted,type='outcome') + ggtitle('Actual outcome of predicted best experiments')
  }
  else
  {
    g_o_p = ggplot() + ggtitle('Actual outcome of predicted best experiments')
  }
  
  g_e = ggarrange(g_e_p, g_e_a, 
                  nrow=1,ncol=2,
                  common.legend = TRUE,
                  align='hv',
                  #labels=c('a', 'b'),
                  legend='bottom')
  
  g_o = ggarrange(g_o_p, g_o_a, 
                  nrow=1,ncol=2,
                  common.legend = TRUE,
                  align='hv',
                  #labels=c('c', 'd'),
                  legend='bottom')
  
  g_final = ggarrange(g_e, g_o, 
                      nrow=2,ncol=1,align='hv')
  cat('+\n')
  
  return(g_final)
}


get_best_case_removal_for_datasets <- function(predictions_best_all)
{
  w = lapply(predictions_best_all, function(dataset_this) { 
    sapply(dataset_this$result_best, function(species_this) { 
      sapply(species_this, function(rep_this) { 
        rep_this$stats$Specificity}) })  } )
  
  result = lapply(1:length(w), function(i) {
    indices = which(w[[i]]==max(w[[i]],na.rm=TRUE),arr.ind=TRUE) %>% head(1) # break ties randomly
    return(list(predictions=predictions_best_all[[i]]$result_best[[ indices[2] ]][[ indices[1] ]],
                index_species=indices[2],
                index_replicate=indices[1]))
  })
  
  return(result)
}

best_case_removals_89 = get_best_case_removal_for_datasets(predictions_best_removal_all_89)
best_case_removals_264 = get_best_case_removal_for_datasets(predictions_best_removal_all_264)

plots_best_removal_89 = lapply(best_case_removals_89, plot_best_experiments_new)
g_best_removal_89 = ggarrange(plotlist=plots_best_removal_89, 
                              nrow=2,ncol=2, labels=sprintf("(%s) %s - species %s, replicate %d", 
                                                                                   letters[1:4],
                                                                                   names_nice[cases$name], 
                                                                                   sapply(best_case_removals_89, function(x) {letters[x$index_species]}),
                                                                                   sapply(best_case_removals_89, function(x) {x$index_replicate})))
ggsave(g_best_removal_89, file='outputs/figures/g_best_removal_89.png',width=20,height=20)

plots_best_removal_264 = lapply(best_case_removals_264, plot_best_experiments_new)
g_best_removal_264 = ggarrange(plotlist=plots_best_removal_264, 
                              nrow=2,ncol=2, labels=sprintf("(%s) %s - species %s, replicate %d", 
                                                            letters[1:4],
                                                            names_nice[cases$name], 
                                                            sapply(best_case_removals_264, function(x) {letters[x$index_species]}),
                                                            sapply(best_case_removals_264, function(x) {x$index_replicate})))
ggsave(g_best_removal_264, file='outputs/figures/g_best_removal_264.png',width=20,height=20)





get_best_case_other_for_datasets <- function(predictions_best_all)
{
  w = lapply(predictions_best_all, function(dataset_this) { 
      sapply(dataset_this, function(rep_this) { 
        rep_this$stats$Specificity}) } )
  print(w)
  
  result = lapply(1:length(w), function(i) {
    indices = which(w[[i]]==max(w[[i]],na.rm=TRUE),arr.ind=TRUE) %>% head(1) # break ties randomly
    return(list(predictions=predictions_best_all[[i]][[ indices[1] ]],
                index_replicate=indices[1]))
  })
  
  return(result)
}


best_case_shannons_h_89 = get_best_case_other_for_datasets(predictions_best_shannons_H_89)
plots_best_shannons_h_89 = lapply(best_case_shannons_h_89, plot_best_experiments_new)
g_best_shannons_h_89 = ggarrange(plotlist=plots_best_shannons_h_89, 
                               nrow=2,ncol=2, labels=sprintf("(%s) %s - replicate %s", 
                                                             letters[1:4],
                                                             names_nice[cases$name], 
                                                             sapply(best_case_shannons_h_89, function(x) {x$index_replicate})))
ggsave(g_best_shannons_h_89, file='outputs/figures/g_best_shannons_h_89.png',width=20,height=20)

best_case_shannons_h_264 = get_best_case_other_for_datasets(predictions_best_shannons_H_264)
plots_best_shannons_h_264 = lapply(best_case_shannons_h_264, plot_best_experiments_new)
g_best_shannons_h_264 = ggarrange(plotlist=plots_best_shannons_h_264, 
                                  nrow=2,ncol=2, labels=sprintf("(%s) %s - replicate %s", 
                                                                letters[1:4],
                                                                names_nice[cases$name], 
                                                                sapply(best_case_shannons_h_264, function(x) {x$index_replicate})))
ggsave(g_best_shannons_h_264, file='outputs/figures/g_best_shannons_h_264.png',width=20,height=20)



best_case_abundance_89 = get_best_case_other_for_datasets(predictions_best_total_abundance_89)
plots_best_abundance_89 = lapply(best_case_abundance_89, plot_best_experiments_new)
g_best_abundance_89 = ggarrange(plotlist=plots_best_abundance_89, 
                                nrow=2,ncol=2, labels=sprintf("(%s) %s - replicate %s", 
                                                              letters[1:4],
                                                              names_nice[cases$name], 
                                                              sapply(best_case_abundance_89, function(x) {x$index_replicate})))
ggsave(g_best_abundance_89, file='outputs/figures/g_best_abundance_89.png',width=20,height=20)

best_case_abundance_264 = get_best_case_other_for_datasets(predictions_best_total_abundance_264)
plots_best_abundance_264 = lapply(best_case_abundance_264, plot_best_experiments_new)
g_best_abundance_264 = ggarrange(plotlist=plots_best_abundance_264, 
                                 nrow=2,ncol=2, labels=sprintf("(%s) %s - replicate %s", 
                                                               letters[1:4],
                                                               names_nice[cases$name], 
                                                               sapply(best_case_abundance_264, function(x) {x$index_replicate})))
ggsave(g_best_abundance_264, file='outputs/figures/g_best_abundance_264.png',width=20,height=20)





# broken below this line
plots_best_shannons_H_264 = lapply(lapply(predictions_best_shannons_H_264, 
                                         pick_highest_sensitivity_replicate_within_dataset),
                                  plot_best_experiments_new)

plots_best_total_abundance_89 = lapply(lapply(predictions_best_total_abundance_89, 
                                           pick_highest_sensitivity_replicate_within_dataset),
                                    plot_best_experiments_new)

plots_best_total_abundance_264 = lapply(lapply(predictions_best_total_abundance_264, 
                                              pick_highest_sensitivity_replicate_within_dataset),
                                       plot_best_experiments_new)

plots_best_removal_89 = lapply(lapply(predictions_best_removal_all_89, # these are already the 1 species we can best find to remove
                                           pick_highest_sensitivity_replicate_within_dataset),
                                    plot_best_experiments_new)

plots_best_removal_264 = lapply(lapply(predictions_best_removal_all_264, # these are already the 1 species we can best find to remove
                                              pick_highest_sensitivity_replicate_within_dataset),
                                       plot_best_experiments_new)

g_best_shannons_H_89 = ggarrange(plotlist=plots_best_shannons_H_89,nrow=2,ncol=2,
                              labels=sprintf("(%s) %s", letters[1:4],names_nice[cases$name]))
g_best_total_abundance_89 = ggarrange(plotlist=plots_best_total_abundance_89,nrow=2,ncol=2,
                                   labels=sprintf("(%s) %s", letters[1:4],names_nice[cases$name]))
g_best_remove_unwanted_89 = ggarrange(plotlist=plots_best_remove_unwanted_89,nrow=2,ncol=2,
                                   labels=sprintf("(%s) %s - Species %s", letters[1:4],names_nice[cases$name], species_this=c(18, 10, 9, 6)))
ggsave(g_best_shannons_H_89, file='outputs/figures/g_best_shannons_H_89.png',width=20,height=20)
ggsave(g_best_total_abundance_89, file='outputs/figures/g_best_total_abundance_89.png',width=20,height=20)
ggsave(g_best_remove_unwanted_89, file='outputs/figures/g_best_remove_unwanted_264.png',width=20,height=20)

g_best_shannons_H_264 = ggarrange(plotlist=plots_best_shannons_H_264,nrow=2,ncol=2,
                                  labels=sprintf("(%s) %s", letters[1:4],names_nice[cases$name]))
g_best_total_abundance_264 = ggarrange(plotlist=plots_best_total_abundance_264,nrow=2,ncol=2,
                                       labels=sprintf("(%s) %s", letters[1:4],names_nice[cases$name]))
g_best_remove_unwanted_264 = ggarrange(plotlist=plots_best_remove_unwanted_264,nrow=2,ncol=2,
                                       labels=sprintf("(%s) %s - Species %s", letters[1:4],names_nice[cases$name], species_this=c(18, 10, 9, 6)))
ggsave(g_best_shannons_H_264, file='outputs/figures/g_best_shannons_H_264.png',width=20,height=20)
ggsave(g_best_total_abundance_264, file='outputs/figures/g_best_total_abundance_264.png',width=20,height=20)
ggsave(g_best_remove_unwanted_264, file='outputs/figures/g_best_remove_unwanted_264.png',width=20,height=20)
















































# do a PCA of the actual outcomes
best_predictions_pca <- function(predictions_this, title="")
{
  df_this = predictions_this$abundances_ALL
  indices_best_actual = predictions_this$indices_best_ACTUAL
  indices_best_predicted = predictions_this$indices_best_PREDICTED
  rows_na = which(is.na(rowSums(df_this)))
  if(length(rows_na)>0)
  {
    df_this = df_this[-rows_na,]
    indices_best_actual = indices_best_actual[-rows_na]
    indices_best_predicted = indices_best_predicted[-rows_na]
  }

  indices_category = factor(paste(indices_best_actual, indices_best_predicted), levels=c("FALSE FALSE","FALSE TRUE", "TRUE FALSE", "TRUE TRUE"), labels=c('true negative','false positive','false negative','true positive'))
  
  print(table(indices_category))
  
  pc_this = prcomp((df_this)^(1/4),scale=TRUE,center=TRUE) # sqrt + centering/scaling
  
  #pc_this$x = pc_this$x[order(indices_category),]
  #indices_category = indices_category[order(indices_category)]
  
  pc_this$rotation = data.frame(pc_this$rotation)
  pc_this$rotation$var = gsub(".star","",row.names(pc_this$rotation), fixed=TRUE)
  
  bin_width_this = max(abs(as.numeric(pc_this$x[,1:2])))
  axis_length_this = max(abs(as.numeric(as.matrix(pc_this$rotation[,1:2]))))
  
  g_biplot = ggplot(data.frame(pc_this$x,indices_category), aes(x=PC1,
                                                                y=PC2)
                                                                ) +
    geom_hex(binwidth=bin_width_this / 5) +
    # geom_point() +
    # geom_jitter(width=0.01,height=0.01) +
    theme_bw() +
    facet_wrap(~indices_category,nrow=1,ncol=4,drop=FALSE) +
    scale_fill_viridis_c(option='plasma',begin=0.1,end=0.9,name='# of cases') +
    geom_segment(data=data.frame(pc_this$rotation),
                 aes(x=0,y=0,xend=PC1*bin_width_this/axis_length_this,yend=PC2*bin_width_this/axis_length_this),
                 arrow=arrow(length = unit(0.05, "inches")),color='darkgray') +
    # geom_text(data=pc_this$rotation,
    #           aes(x=1.1*PC1*bin_width_this/axis_length_this,y=1.1*PC2*bin_width_this/axis_length_this,label=var)) +
    coord_fixed(ratio=1) +
    ggtitle(title)
  
  # g_biplot = ggbiplot(pc_this, 
  #                     groups=indices_category,
  #                     alpha=0.75) +
  #   theme_bw() +
  #   scale_color_manual(values=c('lightblue','orange','red','blue'),drop=FALSE)
  return(g_biplot)
}

best_predictions_pca_all <- function(predictions_raw, is_removal=FALSE)
{
  names_nice_this = names_nice[cases$name]
  g_list = lapply(1:length(predictions_raw), function(i)
  {
    predictions_this = predictions_raw[[i]]$predictions
    
    if (is_removal==FALSE)
    {
      extra_string = ''
    }
    else
    {
      extra_string = sprintf('species %s, ',letters[predictions_raw[[i]]$index_species])
    }
    g = best_predictions_pca(predictions_this, title=sprintf("%s - %sreplicate %d",names_nice[cases$name[i]], extra_string, predictions_raw[[i]]$index_replicate))
  })
  g_final = ggarrange(plotlist = g_list,
                      #align='hv',
                      hjust=0,
                      ncol=1) + ggtitle(title)
  return(g_final)
}

g_best_hexbin_shannons_H_89 = best_predictions_pca_all(best_case_shannons_h_89)
g_best_hexbin_shannons_H_264 = best_predictions_pca_all(best_case_shannons_h_264)
ggsave(ggarrange(g_best_hexbin_shannons_H_89, g_best_hexbin_shannons_H_264,labels='auto'), file='outputs/figures/g_best_hexbin_shannons_H.pdf',width=16,height=10)
ggsave(ggarrange(g_best_hexbin_shannons_H_89, g_best_hexbin_shannons_H_264,labels='auto'), file='outputs/figures/g_best_hexbin_shannons_H.png',width=16,height=10)


g_best_hexbin_abundance_89 = best_predictions_pca_all(best_case_abundance_89)
g_best_hexbin_abundance_264 = best_predictions_pca_all(best_case_abundance_264)
ggsave(ggarrange(g_best_hexbin_abundance_89, g_best_hexbin_abundance_264,labels='auto'), file='outputs/figures/g_best_hexbin_abundance.pdf',width=16,height=10)
ggsave(ggarrange(g_best_hexbin_abundance_89, g_best_hexbin_abundance_264,labels='auto'), file='outputs/figures/g_best_hexbin_abundance.png',width=16,height=10)

g_best_hexbin_removal_89 = best_predictions_pca_all(best_case_removals_89, is_removal = TRUE)
g_best_hexbin_removal_264 = best_predictions_pca_all(best_case_removals_264, is_removal = TRUE)
ggsave(ggarrange(g_best_hexbin_removal_89, g_best_hexbin_removal_264,labels='auto'), file='outputs/figures/g_best_hexbin_removal.pdf',width=16,height=10)
ggsave(ggarrange(g_best_hexbin_removal_89, g_best_hexbin_removal_264,labels='auto'), file='outputs/figures/g_best_hexbin_removal.png',width=16,height=10)

