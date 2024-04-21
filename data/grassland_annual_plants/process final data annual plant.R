library(dplyr)

df_raw = read.csv('assemblages_annual_plant_18.csv') %>%
  select(-stable, -feasible, -richness)

df_names = read.csv('speciesvitalrates.csv')

names_all = df_names %>% 
  pull(species)

names(df_raw)[1:length(names_all)] = paste(names_all,"action",sep=".")
names(df_raw)[which(grepl("star$",names(df_raw)))] = paste(names_all,"outcome",sep=".")
      
df_final = df_raw %>%
  mutate(environment.initial = 0) %>%
  select(contains("action"),contains("initial"),contains("outcome"))

write.csv(df_final, file='data_grassland_annual_plants.csv', row.names = FALSE)
