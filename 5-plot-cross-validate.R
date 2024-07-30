library(ggplot2)
library(dplyr)
library(stringr)


# Setup output directory
try(dir.create(file.path(getwd(), 'outputs/figures'), recursive = TRUE))
try(dir.create(file.path(getwd(), 'outputs/statistical'), recursive = TRUE))
directory_string = file.path(getwd(), 'outputs/statistical')


results_all = read.csv('outputs/statistical/cross_validation.csv') %>%
  mutate(name_nice=str_to_sentence(gsub("_", " ", results_all$name))) %>%
  mutate(name_nice = gsub("Grassland annual","Grassland annual\n",name_nice)) %>%
  mutate(test_mode_nice = gsub(" environments","",gsub("_"," ",gsub("mae_test_","",test_mode))))



g_cross_validation = ggplot(results_all, aes(x=num_train_cases,y=value_scaled,color=test_mode_nice)) + 
  geom_point(alpha=0.5) +
  #scale_color_viridis_d() +
  facet_grid(name_nice~num_train_env) +
  theme_bw() + 
  xlab('Number of training cases') +
  ylab('Scaled error') + 
  geom_smooth(method='lm') +
  scale_x_log10() +
  #scale_y_log10() +
  scale_y_log10(breaks=c(0.005,0.01,0.05,0.2),limits=c(0.005,0.2)) +
  scale_color_brewer(name='Test environments',palette='Set1')

ggsave(g_cross_validation, file='outputs/figures/g_cross_validation.png',width=8,height=4.5)
