library(dplyr)
library(data.table)

data_fly = read.csv('FlygutCFUsData.csv')

data_initial = data_fly %>% 
  select(LP:AO)
names(data_initial) = letters[1:ncol(data_initial)]

data_final = data_fly %>%
  select(LP.abund..CFU. : AO.abund..CFU.)
names(data_final) = paste(letters[1:ncol(data_final)],"star",sep=".")

richness = apply(data_final, 1, function(x) {sum(x>0)})
stable = rep(NA, nrow(data_final))
feasible = rep(TRUE, nrow(data_final))

data_export = data.frame(data_initial, stable=stable, feasible=feasible, richness=richness, data_final)

write.csv(data_export, file='data_fly.csv',row.names=F)
