library(dplyr)
library(data.table)

data_fly = read.csv('FlygutCFUsData.csv')

data_action = data_fly %>% 
  select(LP:AO)
names(data_action) = paste(names(data_action), "action", sep=".")

data_outcome = data_fly %>%
  select(LP.abund..CFU. : AO.abund..CFU.)
names(data_outcome) = gsub(".abund..CFU.","",names(data_outcome))
names(data_outcome) = paste(names(data_outcome), "outcome", sep=".")

data_export = data.frame(data_action, environment.initial=0, data_outcome)

write.csv(data_export, file='data_fly_gut.csv',row.names=F)
