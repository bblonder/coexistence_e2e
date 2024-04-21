library(ggplot2)
library(tidyr)
library(dplyr)


# data are from https://github.com/pennekampster/Code_and_data_OverallEcosystemStability

data = read.csv('species_biomass_BEEP_OES.csv')


# 
# data_ss = data[1:1000,]
# 
# ggplot(data_ss, aes(x=day, y=species_biomass,color=predicted_species)) + 
#   geom_line() +
#   facet_wrap(~paste(combination, temperature, replicate))

data_wide = data %>%
  group_by(combination, temperature, replicate, predicted_species) %>%
  # pick the last day of the experiment
  slice_max(day,n=1) %>%
  select(-X, -microcosmID, -day, -total_biomass, -richness) %>%
  # put the species biomass as the outcome
  pivot_wider(names_from='predicted_species', values_from='species_biomass',names_glue = "{predicted_species}.outcome")

# for each combination, fill in the presence/absence matrix
combinations = data.frame(combination=unique(data_wide$combination))

names_split = strsplit(combinations$combination,", ")
index_longest = which.max(sapply(names_split, length))
names_all = names_split[[index_longest]]
combinations = cbind(combinations,data.frame(matrix(data=0,nrow=nrow(combinations),ncol=length(names_all),dimnames=list(NULL, names_split[[index_longest]]))))

for (i in 1:length(names_split))
{
  combinations[i, names_split[[i]]] = 1
}

# join the presence/absence to the experimental data
data_wide_processed = data_wide %>% 
  left_join(combinations, by='combination') %>%
  ungroup %>%
  select(-combination, -replicate) %>%
  # replace all NAs with 0s
  replace(is.na(.), 0) %>%
  # rearrange
  select(all_of(names_all), temperature.initial = temperature, all_of(paste(names_all,"outcome",sep=".")),)
names(data_wide_processed)[1:length(names_all)] = paste(names(data_wide_processed)[1:length(names_all)], "action",sep=".")

write.csv(data_wide_processed, file='data_ciliates.csv', row.names=FALSE)
