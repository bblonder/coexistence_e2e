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
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")



colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


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
                  num_species_dataset='Dataset, regional pool # species',
                  type='Dataset type',
                  num_losses_mean_train='Experiment, ∆ species richness',
                  abundance_final_skewness_mean_train='Outcome, abundance skewness',
                  nice_name='Dataset')

experimental_design_nice = c(prior='Singlets + 1-dropouts, then mixed',
                             mixed='Mixed',
                             `low-2`='Doublets only',
                             `low-3`='Doublets + triplets only',
                             `high-1`='1-dropouts only',
                             `high-2`='1-dropouts + 2-dropouts only')

methods_nice = c(naive=' Naïve (mean abundance) ',
                             rf='Random forest on experiments',
                             glv='GLV predictions',
                             glv_rf='Random forest on GLV residuals',
                             glv_rf_full='Random forest on experiments + GLV predictions',
                             sequential_rf='Random forest, sequential')

# add NA removal counts
source('utils/quantile_trim.R')
fns = c(`cedar_creek_plants`='data/cedar_creek/cedar_creek_2018.csv', 
        `sortie-nd_plants`='data/sortie/data_sortie.csv',
        `human_gut`='data/glv/assemblages_H_12.csv',
        `mouse_gut`='data/glv/assemblages_M_11.csv',
        `glv_simulated`='data/glv/assemblages_glv_16.csv',
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
  mutate(experimental_design_nice = factor(experimental_design_nice[experimental_design])) %>%
  mutate(method_nice = factor(methods_nice[method])) %>%  
  
  mutate(empirical=name %in% c("cedar_creek_plants","fly_gut","soil_bacteria")) %>%
  mutate(num_na = num_nas[name]) %>%
  mutate(num_states = row_counts[name]) %>%
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
         num_states, 
         abundance_q95_dataset,
         num_na, 
         empirical) %>%
  unique %>%
  arrange(name_nice)
write.csv(df_all_stats,'outputs/figures/table_dataset_stats.csv',row.names=F)











### PERFORMANCE OVERALL

#,levels=c('prior','mixed','low-2','low-3','high-1','high-2')





# ggplot(df_all_for_regression %>%
#          filter(num_train==30), 
#        aes(x=experimental_design,
#            y=abundance_mae_mean_test_scaled_clipped,
#            color=name)) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_wrap(~method) +
#   ylim(0,1)


# # good plot
# ggplot(df_all %>%
#          filter(experimental_design=='mixed'), 
#        aes(x=method_nice,
#            y=abundance_mae_mean_test_scaled_clipped,
#            color=factor(num_train))) +
#   geom_boxplot() +
#   theme_bw() +
#   facet_wrap(~name_nice) +
#   ylim(0,1) +
#   #scale_x_discrete(drop=FALSE) +
#   scale_color_viridis_d() +
#   coord_flip()


ggplot(df_all %>%
         filter(experimental_design=='mixed'), 
       aes(x=num_train,
           y=abundance_mae_mean_test_scaled_clipped,
           color=name_nice,fill=name_nice)) +
  geom_point(alpha=0.25) +
  theme_bw() +
  facet_wrap(~method_nice) +
  ylim(0,10) +
  #scale_x_discrete(drop=FALSE) +
  scale_color_brewer(palette='Set2') +
  scale_fill_brewer(palette='Set2') +
  geom_smooth(method='lm') +
  scale_x_log10() +
  scale_y_sqrt(breaks=c(0,0.1,1,5,10))




# good plot
ggplot(df_all_for_regression %>%
         filter(num_train==30), 
       aes(x=method,
           y=abundance_mae_mean_test_scaled_clipped,
           color=experimental_design_nice)) +
  geom_boxplot() +
  theme_bw() +
  ylim(0,1) +
  facet_wrap(~name) +
  #scale_x_discrete(drop=FALSE) +
  scale_color_brewer(palette='Set1')




# df_all_for_regression %>%
#   filter(num_train==30) %>%
#   filter(name=='Annual plant') %>%
#   filter(experimental_design=='mixed') %>%
#   filter(method=='glv_rf')
# 








# fit models
m_richness_scaled_mae_test = lmer(richness_scaled_mae_test~method*log10(num_train)*experimental_design*name + (1|name),
                                  data=df_all_for_regression)
m_composition_balanced_accuracy_mean_test = lmer(oneminuscomposition_balanced_accuracy_mean_test~method*log10(num_train)*experimental_design*name + (1|name),
                                                 data=df_all_for_regression)
m_abundance_scaled_mae_test = lmer(abundance_scaled_mae_test~method*log10(num_train)*experimental_design*name + (1|name),
                                   data=df_all_for_regression %>% filter(method!="sequential_rf"))

r.squaredGLMM(m_richness_scaled_mae_test)
r.squaredGLMM(m_composition_balanced_accuracy_mean_test)
r.squaredGLMM(m_abundance_scaled_mae_test)


# q1a result
plot_visreg_by_method <- function(model, ylab, title)
{
  g_method = visreg(model, xvar='name',by='method',gg=TRUE,
                    cond=list(num_train=100,experimental_design='mixed'), 
                    overlay=TRUE, band=FALSE,
                    points=list(size=0.1)
  ) +
    theme_bw() +
    xlab('Dataset') +
    ylab(ylab) +
    labs(color='Method') +
    scale_color_manual(values=colorBlindBlack8) +
    ggtitle(title) +
    ylim(0,1)
}

g_visreg_richness_method = plot_visreg_by_method(m_richness_scaled_mae_test, "Scaled mean absolute error",title='Richness')
g_visreg_composition_method = plot_visreg_by_method(m_composition_balanced_accuracy_mean_test, "1 - balanced accuracy",title='Composition')
g_visreg_abundance_method = plot_visreg_by_method(m_abundance_scaled_mae_test, "Scaled mean absolute error",title='Abundance')
g_visreg_by_method = ggarrange(g_visreg_richness_method, 
                               g_visreg_composition_method,
                               g_visreg_abundance_method,
                               nrow=3,ncol=1,
                               labels='auto',
                               common.legend = TRUE,
                               align='hv',
                               legend='bottom')
ggsave(g_visreg_by_method, 
       file='outputs/figures/g_visreg_by_method.png',
       width=7,height=8)

plot_visreg_by_numtrain <- function(model, ylab, title)
{
  g_numtrain = visreg(model, xvar='num_train',by='name',gg=TRUE,
                      cond=list(method='love',experimental_design='mixed'),
                      overlay=TRUE, band=FALSE,
                      points=list(size=0.1)) +
    theme_bw()+
    scale_x_log10() +
    xlab('Number of training cases') +
    ylab(ylab) +
    labs(color='Dataset') +
    scale_color_manual(values=colorBlindBlack8) +
    ggtitle(title) +
    ylim(0,1)
}

g_visreg_richness_numtrain = plot_visreg_by_numtrain(m_richness_scaled_mae_test, "Scaled mean absolute error", "Richness")
g_visreg_composition_numtrain = plot_visreg_by_numtrain(m_composition_balanced_accuracy_mean_test, "1 - balanced accuracy", "Composition")
g_visreg_abundance_numtrain = plot_visreg_by_numtrain(m_abundance_scaled_mae_test, "Scaled mean absolute error", "Abundance")
g_visreg_by_numtrain = ggarrange(g_visreg_richness_numtrain,
                                 g_visreg_composition_numtrain,
                                 g_visreg_abundance_numtrain,
                                 nrow=3,ncol=1,
                                 labels='auto',
                                 common.legend = TRUE,
                                 align='hv',
                                 legend='bottom')
ggsave(g_visreg_by_numtrain, 
       file='outputs/figures/g_visreg_by_numtrain.png',
       width=7,height=8)



# now look at experimental design
plot_visreg_by_experimental_design <- function(model, ylab, title)
{
  g_experimental_design = visreg(model, xvar='name',by='experimental_design',gg=TRUE,
                                 cond=list(method='love',num_train=100),
                                 overlay=TRUE, band=FALSE,
                                 points=list(size=0.1)) +
    theme_bw() +
    xlab('Dataset') +
    ylab(ylab) +
    labs(color='Experimental design') +
    scale_color_manual(values=colorBlindBlack8) +
    ggtitle(title) +
    ylim(0,1)
}
g_visreg_richness_experimental_design = plot_visreg_by_experimental_design(m_richness_scaled_mae_test, "Scaled mean absolute error", "Richness")
g_visreg_composition_experimental_design = plot_visreg_by_experimental_design(m_composition_balanced_accuracy_mean_test, "1 - balanced accuracy", "Composition")
g_visreg_abundance_experimental_design = plot_visreg_by_experimental_design(m_abundance_scaled_mae_test, "Scaled mean absolute error", "Abundance")
g_visreg_by_experimental_design = ggarrange(g_visreg_richness_experimental_design, 
                                            g_visreg_composition_experimental_design,
                                            g_visreg_abundance_experimental_design,
                                            nrow=3,ncol=1,
                                            labels='auto',
                                            common.legend = TRUE,
                                            align='hv',
                                            legend='bottom')
ggsave(g_visreg_by_experimental_design, file='outputs/figures/g_visreg_by_experimental_design.png',width=7,height=8)





### explanations
df_all_for_regression_dataset = df_all_for_regression %>%
  filter(num_train==89 & method=='love') %>%
  arrange(num_train, experimental_design, rep) %>%
  mutate(type = factor(empirical,levels=c(FALSE,TRUE),labels = c('Simulated','Empirical')))

xvars_all = c("type",
              "num_species_dataset", 
              "num_losses_mean_train",
              "abundance_final_skewness_mean_train")

m_richness_scaled_mae_test_dataset = lmer(richness_scaled_mae_test ~ (type + num_species_dataset + num_losses_mean_train + abundance_final_skewness_mean_train) + (1|nice_name),
                                          data=df_all_for_regression_dataset)

m_composition_balanced_accuracy_mean_test_dataset = lmer(oneminuscomposition_balanced_accuracy_mean_test ~ (type + num_species_dataset + num_losses_mean_train + abundance_final_skewness_mean_train) + (1|nice_name),
                                                         data=df_all_for_regression_dataset)

m_abundance_scaled_mae_test_dataset = lmer(abundance_scaled_mae_test ~ (type + num_species_dataset + num_losses_mean_train + abundance_final_skewness_mean_train) + (1|nice_name),
                                           data=df_all_for_regression_dataset)



plot_visreg_perf_dataset <- function(data=df_all_for_regression_dataset, 
                                     yvar, 
                                     xvars=c('type','num_species_dataset','num_losses_mean_train', 'abundance_final_skewness_mean_train'),
                                     ylab,
                                     ylim,
                                     no_legend=FALSE)
{
  model = lmer(formula(sprintf('%s ~ %s + (num_losses_mean_train|nice_name)', yvar, paste(xvars,collapse=" + "))),
               data=data)
  
  plots_all = vector(mode="list",length=length(xvars))
  for (i in 1:length(xvars))
  {
    plots_all[[i]] <- visreg(model, xvar=xvars[i], 
                             gg=TRUE, overlay=TRUE, band=FALSE, 
                             by='nice_name',
                             line=list(col="black")) +
      xlab(varnames_nice[ xvars[i] ]) +
      theme_bw() +
      ylab(ylab) +
      ylim(ylim[1],ylim[2]) +
      labs(color='Dataset') + 
      scale_color_manual(values=colorBlindBlack8) + theme(plot.margin = margin(t=0.7,0.1,0.1,0.1, "cm"))
  }
  #plots_arranged = ggarrange(plotlist=plots_all, nrow=1,ncol=4,common.legend = TRUE, legend=ifelse(no_legend==TRUE,'none','bottom'),align='hv') +
  #  theme(plot.margin = margin(0.5,0.1,0.1,0.1, "cm"))
  #plots_arranged = annotate_figure(plots_arranged, fig.lab=ylab)
  
  return(list(model=model,plots=plots_all))
}

perf_dataset_richness = plot_visreg_perf_dataset(yvar='richness_scaled_mae_test',
                                                 ylab='Scaled mean absolute error',
                                                 ylim=c(0,0.25),
                                                 no_legend = FALSE)

perf_dataset_composition = plot_visreg_perf_dataset(yvar='oneminuscomposition_balanced_accuracy_mean_test',
                                                    ylab='One minus balanced accuracy',
                                                    ylim=c(0,0.8),
                                                    no_legend = FALSE)

perf_dataset_abundance = plot_visreg_perf_dataset(yvar='abundance_scaled_mae_test',
                                                  ylab='Scaled mean absolute error',
                                                  ylim=c(0,0.4),
                                                  no_legend = FALSE)

g_perf_dataset = ggarrange(plotlist=c(perf_dataset_richness$plots, perf_dataset_composition$plots, perf_dataset_abundance$plots),
                           nrow=3,ncol=4,
                           common.legend = TRUE,
                           align='hv',
                           legend='bottom',
                           labels=c('(a) Richness',rep("",3),'(b) Composition',rep("",3),'(c) Abundance',rep("",3)))

ggsave(g_perf_dataset, file='outputs/figures/g_perf_dataset.png',
       width=11,height=9)

r.squaredGLMM(perf_dataset_richness$model)
r.squaredGLMM(perf_dataset_composition$model)
r.squaredGLMM(perf_dataset_abundance$model)





#### FIND BEST CASES
source('pick_datasets.R')
source('log_seq.R')
possible_num_train = ceiling(log_seq(1e1,1e4,length.out=20))

cases = data.frame(name=names(nice_names[nice_names != "GLV simulated"]), fn_assemblages=c('data/annual_plant/assemblages_annual_plant_18.csv',
                                                                                           'data/cedar_creek/cedar_creek_2018.csv',
                                                                                           'data/fly/data_fly.csv',
                                                                                           #'data/glv/assemblages_glv_16.csv',
                                                                                           'data/glv/assemblages_H_12.csv',
                                                                                           'data/glv/assemblages_M_11.csv',
                                                                                           'data/sortie/data_sortie.csv',
                                                                                           'data/friedman_gore/data_friedman_gore.csv'))

# only do the complete sets
dataset_ids_best_experiments = which(cases$name %in% c('annual_plant','human_gut','mouse_gut','sortie-nd_plants'))


datasets_preds = lapply(1:nrow(cases), function(i) {
  cat(i)
  datasets_preds_this = pick_datasets(name=cases$name[i], 
                                      method='love',
                                      num_train=89,
                                      experimental_design = 'mixed',
                                      response_var = 'abundance')
})
names(datasets_preds) = cases$name


predictions_best <- function(type, datasets_preds, quantile_cutoff, fn_assemblages, options=NULL)
{
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
    
    if (type=='shannons_H')
    {
      shannons_H_PREDICTED = diversity(outcome_abundances_PREDICTED,index='shannon')
      shannons_H_ALL = diversity(outcome_abundances_ALL,index='shannon')
      
      experiments_best_PREDICTED = experiments[shannons_H_PREDICTED >= quantile(shannons_H_ALL, quantile_cutoff, na.rm=T), 1:num_species] %>%
        na.omit
      experiments_best_ACTUAL = assemblages[shannons_H_ALL >= quantile(shannons_H_ALL, quantile_cutoff, na.rm=T), 1:num_species] %>%
        na.omit
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
    
    return(list(name=datasets_preds$name[1],
                type=type,
                num_species=num_species,
                stats=confusion_matrix_stats,
                experiments_ALL=ids_experiments_ALL,
                experiments_best_ACTUAL=ids_experiments_best_ACTUAL,
                experiments_best_PREDICTED=ids_experiments_best_PREDICTED))
  })
  cat('\n')
  return(result)
}

check_unwanted <- function(dataset_name_this)
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
    
    result_this = rbindlist(lapply(predictions_best_remove_unwanted, function(x) { x$stats }))
    result_this$id_unwanted = id_unwanted_this
    result_this$name = dataset_name_this
    return(list(result_this=result_this,predictions_best_remove_unwanted=predictions_best_remove_unwanted))
  })
  
  result_df = rbindlist(lapply(stats_unwanted_this, function(x) { x$result_this }))
  result_best = lapply(stats_unwanted_this, function(x) { x$predictions_best_remove_unwanted })
  return(list(result_df=result_df, result_best=result_best))
}
results_unwanted_all = lapply(cases$name, check_unwanted) # this runs slowly due to the annual plant dataset

# only show the complete datasets
g_unwanted_all = ggplot(lapply(results_unwanted_all[dataset_ids_best_experiments], function(x) {x[[1]]}) %>% 
                          rbindlist %>% 
                          select(name, id_unwanted, Sensitivity, Specificity) %>% 
                          melt(id.vars=c("name","id_unwanted")) %>%
                          mutate(value=ifelse(is.na(value),0,value)) %>%
                          mutate(id_unwanted = letters[id_unwanted]), 
                        aes(x=factor(id_unwanted),y=value,color=variable)) + 
  facet_wrap(~name,scales='free_x',labeller=as_labeller(nice_names),nrow=1) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Remove unwanted species") +
  xlab("Unwanted species") + ylab("Value (train + test)") +
  scale_color_manual(values=colorBlindBlack8,name='Statistic') +
  theme(legend.position='bottom') +
  theme(axis.text=element_text(size=6))




predictions_best_shannons_H = lapply(1:length(datasets_preds),
                                     FUN=function(i)
                                     {
                                       predictions_best(datasets_preds=datasets_preds[[i]],
                                                        type='shannons_H',
                                                        quantile_cutoff=0.95,
                                                        fn_assemblages=cases$fn_assemblages[i])
                                     })
predictions_best_total_abundance = lapply(1:length(datasets_preds),
                                          FUN=function(i)
                                          {
                                            predictions_best(datasets_preds=datasets_preds[[i]],
                                                             type='total_abundance',
                                                             quantile_cutoff=0.95,
                                                             fn_assemblages=cases$fn_assemblages[i])
                                          })






# plot performance
best_experiments_classification_statistics <- function(predictions_best, title)
{
  df_best_experiments_classification_statistics = lapply(predictions_best, 
                                                         FUN=function(predictions_this)
                                                         {
                                                           dataset_name = predictions_this[[1]]$name
                                                           
                                                           stats_this = rbindlist(lapply(predictions_this, function(x) { x$stats }))
                                                           stats_this$name = dataset_name
                                                           return(stats_this)
                                                         }) %>% 
    rbindlist %>% 
    select(name,Sensitivity,Specificity) %>% 
    melt(id.vars='name')
  
  df_best_experiments_classification_statistics$value[is.na(df_best_experiments_classification_statistics$value)] = 0
  
  ggplot(df_best_experiments_classification_statistics %>% 
           mutate(nice_name=nice_names[name]), 
         aes(x=variable,y=value,color=nice_name)) +
    geom_boxplot() +
    ylim(0,1) +
    theme_bw() +
    xlab("Statistic") +
    ylab("Value (train + test)") +
    scale_color_manual(values=colorBlindBlack8,name='Dataset') +
    ggtitle(title)
}


g_best_experiments_classification_statistics_shannons_H = best_experiments_classification_statistics(predictions_best_shannons_H[dataset_ids_best_experiments], 
                                                                                                     ">95% quantile Shannon's H")
g_best_experiments_classification_statistics_total_abundance = best_experiments_classification_statistics(predictions_best_total_abundance[dataset_ids_best_experiments], 
                                                                                                          ">95% quantile total abundance")
ggsave(ggarrange(g_unwanted_all,
                 ggarrange(g_best_experiments_classification_statistics_shannons_H, 
                           g_best_experiments_classification_statistics_total_abundance,
                           nrow=1,ncol=2,labels=c('b','c'),align='hv',common.legend = TRUE, legend='bottom'),
                 nrow=2,ncol=1,labels=c('a','')),
       file='outputs/figures/g_best_experiments_classification_statistics.png',
       width=8,height=6)



prep_best_experiments_matrix <- function(fn_assemblages, predictions_this)
{
  num_species = predictions_this[[1]]$num_species
  
  #experiments_all = read.csv(fn_assemblages)[1:num_species]
  #ids_experiments_all = apply(experiments_all, 1, paste, collapse="*")
  df_all_reps = rbindlist(lapply(1:length(predictions_this), function(i)
  {
    cat('.')
    ids_predicted_rep = predictions_this[[i]]$experiments_best_PREDICTED
    ids_actual_rep = predictions_this[[i]]$experiments_best_ACTUAL
    ids_all_rep = predictions_this[[i]]$experiments_ALL
    
    df_rep = data.frame(rep=i, 
                        id=ids_all_rep,
                        predicted = ids_all_rep %in% ids_predicted_rep,
                        actual = ids_all_rep %in% ids_actual_rep)
    return(df_rep)
  })) %>%
    mutate(true.positive=(predicted==TRUE & actual==TRUE),
           true.negative=(predicted==FALSE & actual==FALSE),
           false.positive=(predicted==TRUE & actual==FALSE),
           false.negative=(predicted==FALSE & actual==TRUE))
  df_all_reps$classification = apply(df_all_reps %>% select(true.positive:false.negative), 1, function(x) {gsub("\\.", " ", names(x)[which(x==TRUE)])})
  df_all_reps$classification = factor(df_all_reps$classification, levels=c('false negative','false positive','true negative','true positive'))
  df_all_reps$name = predictions_this[[1]]$name
  df_all_reps$type = predictions_this[[1]]$type
  df_all_reps$extra_text = ""
  
  return(df_all_reps)
}


best_experiments_matrix_shannons_H = lapply(1:nrow(cases), FUN=function(i) {
  cat('+')
  prep_best_experiments_matrix(fn_assemblages = cases$fn_assemblages[i], 
                               predictions_this = predictions_best_shannons_H[[i]])
})

best_experiments_matrix_total_abundance = lapply(1:nrow(cases), FUN=function(i) {
  cat('+')
  prep_best_experiments_matrix(fn_assemblages = cases$fn_assemblages[i], 
                               predictions_this = predictions_best_total_abundance[[i]])
})

best_rows_unwanted = lapply(results_unwanted_all, function(x) {x[[1]] %>% 
    select(name, id_unwanted, Specificity)}) %>% 
  rbindlist %>% 
  group_by(name, id_unwanted) %>% 
  mutate(Specificity = ifelse(is.na(Specificity),0,Specificity)) %>%
  summarize(Specificity.mean=mean(Specificity, na.rm=T)) %>%
  filter(Specificity.mean==max(Specificity.mean,na.rm=T)) %>%
  filter(row_number()==1)
# reorder by case
best_rows_unwanted = best_rows_unwanted[match(best_rows_unwanted$name, cases$name),]

# pick bets experiments as those with maximum specificity
predictions_best_remove_unwanted = lapply(1:nrow(best_rows_unwanted), function(dataset_index) {
  results_unwanted_all[[ dataset_index ]][[ 2 ]][[ best_rows_unwanted$id_unwanted[dataset_index] ]]
})

best_experiments_matrix_remove_unwanted = lapply(1:nrow(cases), FUN=function(i) {
  cat('+')
  prep_best_experiments_matrix(fn_assemblages = cases$fn_assemblages[i],
                               predictions_this = predictions_best_remove_unwanted[[i]])
})




# make plots of best experiments


plot_best_experiments_error_rates <- function(df, extra_text="")
{
  df_short = df %>% 
    group_split(rep) %>% 
    #lapply(function(x) {x$actual}) %>%
    lapply(function(x) {x$classification}) %>% 
    as.data.frame
  names(df_short) = 1:ncol(df_short)
  df_long = reshape2::melt(as.matrix(df_short)) %>%
    mutate(Var2=as.character(Var2)) %>%
    #mutate(value=factor(value))
    mutate(value=factor(value,levels=levels(df$classification))) %>%
    mutate(Var2=factor(Var2,levels=1:10,labels=paste(1:10),ordered=TRUE))
  
  g = ggplot(df_long,aes(x=Var1,y=Var2,fill=value)) + 
    geom_tile() + 
    scale_fill_manual(values=c('red','orange3','black','blue'),name='Classification',drop=FALSE) + 
    theme_bw() +
    xlab('Experiment') +
    ylab('Replicate') +
    ggtitle(paste(df$name[1], df$type[1])) +
    scale_x_continuous(expand=c(0,0),breaks=range(df_long$Var1)) +
    ggtitle(paste(nice_names[df$name[1]], extra_text)) + 
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text.x=element_text(size=3))
  
  return(g)
}

plots_best_experiments_error_rates_shannons_H = lapply(dataset_ids_best_experiments, function(i) {
  plot_best_experiments_error_rates(best_experiments_matrix_shannons_H[[i]], 
                                    extra_text = " - find >95% quantile Shannon's H")
})

plots_best_experiments_error_rates_total_abundance = lapply(dataset_ids_best_experiments, function(i) {
  plot_best_experiments_error_rates(best_experiments_matrix_total_abundance[[i]], 
                                    extra_text = " - find >95% quantile total abundance")
})

plots_best_experiments_error_rates_remove_unwanted = lapply(dataset_ids_best_experiments, function(i) {
  plot_best_experiments_error_rates(best_experiments_matrix_remove_unwanted[[i]], 
                                    extra_text = sprintf(" - remove %s", letters[ best_rows_unwanted$id_unwanted[i] ]))
})

g_best_experiments_error_rates_all = ggarrange(plotlist=c(plots_best_experiments_error_rates_remove_unwanted,
                                                          plots_best_experiments_error_rates_shannons_H,
                                                          plots_best_experiments_error_rates_total_abundance
),
nrow=3,ncol=4,
common.legend = TRUE,
labels=c('(a)',rep('',3),'(b)',rep('',3),'(c)',rep('',3)),
legend='bottom',
align='hv')
ggsave(g_best_experiments_error_rates_all,
       file='outputs/figures/g_best_experiments_error_rates_all.png',
       width=12,height=8,dpi=1200)





plot_best_experiments_presence_absence <- function(fn_assemblages, df, extra_text="")
{
  cat("*")
  
  experiments_all = read.csv(fn_assemblages) 
  num_species = ncol(experiments_all %>% select(contains("star")))
  experiments_all = experiments_all[,1:num_species]
  ids_experiments_all = apply(experiments_all, 1, paste, collapse="*")
  
  experiments_to_plot = experiments_all
  # convert to character format
  for (j in 1:ncol(experiments_to_plot))
  {
    experiments_to_plot[,j] = as.character(experiments_to_plot[,j])
  }
  
  # pick the best output
  df_count_pos = df %>% 
    group_by(rep) %>% 
    summarize(tp=sum(classification=='true positive')) %>%
    ungroup
  
  id_best = df_count_pos$rep[which.max(df_count_pos$tp)]
  #print(df_count_pos)
  
  df_this = df %>% filter(rep==id_best)
  
  # make sure the ordering of rows matches
  stopifnot(all((df %>% filter(rep==id_best) %>% pull(id)) == ids_experiments_all))
  
  experiments_to_plot_final = matrix(rep(as.character(df_this$classification),ncol(experiments_to_plot)),
                                     nrow=nrow(experiments_to_plot),ncol=ncol(experiments_to_plot))
  experiments_to_plot_final[experiments_to_plot=="0"] = NA
  
  experiments_to_plot_long = reshape2::melt(as.matrix(experiments_to_plot_final)) %>%
    mutate(value=factor(value,levels=levels(df$classification))) %>%
    mutate(Var2=factor(Var2,levels=1:num_species,labels=letters[1:num_species]))
  
  g = ggplot(experiments_to_plot_long,aes(x=Var1,y=Var2,fill=value)) + 
    geom_tile() + 
    scale_fill_manual(values=c('red','orange3','black','blue'),name='Classification',drop=FALSE) +
    theme_bw() +
    xlab('Experiment') +
    ylab('Species') +
    scale_x_continuous(expand=c(0,0),breaks=range(experiments_to_plot_long$Var1)) +
    ggtitle(paste(nice_names[df$name[1]], "-", "replicate", id_best, "-", extra_text)) + 
    theme(plot.title = element_text(size = 8)) +
    theme(axis.text.x=element_text(size=3))
  
  cat('\n')
  
  return(g)
}

plots_best_experiments_presence_absence_shannons_H = lapply(dataset_ids_best_experiments,
                                                            FUN=function(i)
                                                            {
                                                              plot_best_experiments_presence_absence(
                                                                fn_assemblages=cases$fn_assemblages[i],
                                                                df=best_experiments_matrix_shannons_H[[i]],
                                                                extra_text = "\nfind >95% quantile Shannon's H")
                                                            })

plots_best_experiments_presence_absence_total_abundance = lapply(dataset_ids_best_experiments,
                                                                 FUN=function(i)
                                                                 {
                                                                   plot_best_experiments_presence_absence(
                                                                     fn_assemblages=cases$fn_assemblages[i],
                                                                     df=best_experiments_matrix_total_abundance[[i]],
                                                                     extra_text = "\nfind >95% quantile total abundance")
                                                                 })

plots_best_experiments_presence_absence_remove_unwanted = lapply(dataset_ids_best_experiments,
                                                                 FUN=function(i)
                                                                 {
                                                                   plot_best_experiments_presence_absence(
                                                                     fn_assemblages=cases$fn_assemblages[i],
                                                                     df=best_experiments_matrix_remove_unwanted[[i]],
                                                                     extra_text = sprintf("\nremove %s", letters[ best_rows_unwanted$id_unwanted[i] ]))
                                                                 })


g_best_experiments_presence_absence_all = ggarrange(plotlist=c(plots_best_experiments_presence_absence_remove_unwanted,
                                                               plots_best_experiments_presence_absence_shannons_H,
                                                               plots_best_experiments_presence_absence_total_abundance
),
nrow=3,ncol=4,
common.legend = TRUE,

labels=c('(a)',rep('',3),'(b)',rep('',3),'(c)',rep('',3)),
legend='bottom',
align='hv'
)

ggsave(g_best_experiments_presence_absence_all,
       file='outputs/figures/g_best_experiments_presence_absence_all.png',
       width=12,height=8,dpi=1200)



















### END NEW CODE SECTION


source('pick_datasets.R')

predictions_abundance_all_for_scatter = rbindlist(lapply(1:nrow(cases), function(j) {
  cat('.')
  z = pick_datasets(name=cases$name[j], 
                    method='love',
                    num_train=89,
                    experimental_design = 'mixed',
                    response_var = 'abundance')
  
  z_all = rbindlist(lapply(1:length(z$files_pred_test), function(i)
  {
    df_this = cbind(reshape2::melt(read.csv(z$files_pred_test[i])) %>%
                      rename(abundance_pred=value),
                    reshape2::melt(read.csv(z$files_obs_test[i])) %>%
                      rename(abundance_obs=value) %>%
                      select(-variable)) %>%
      mutate(rep=i) %>%
      mutate(variable=gsub("\\.star","",variable))
    return(df_this)
  }))
  
  z_all$name = cases$name[j]
  
  return(z_all)
}))

# reduce # of cases to make the plotting work
predictions_abundance_all_for_scatter_small = predictions_abundance_all_for_scatter %>%
  group_by(name, rep, variable) %>%
  sample_n(size=min(100,n())) %>%
  mutate(rep.factor=factor(rep))

g_abundance_scatter = ggplot(predictions_abundance_all_for_scatter_small, aes(x=abundance_pred,y=abundance_obs,color=variable,shape=rep.factor)) + 
  geom_line(stat="smooth",method = "lm", se=FALSE,alpha=0.5,size=0.75) +
  geom_point() +
  theme_bw() +
  geom_abline(slope=1,size=2) +
  scale_shape_manual(values=1:nlevels(predictions_abundance_all_for_scatter_small$rep.factor)) + 
  scale_color_hue(l=50,h=c(10,350)) +
  facet_wrap(~name,scales='free',labeller=as_labeller(nice_names),nrow=2) +
  labs(color='Species',shape='Replicate') +
  xlab('Abundance (predicted)') +
  ylab('Abundance (observed)')

ggsave(g_abundance_scatter, file='outputs/figures/g_abundance_scatter.png',width=12,height=8)


