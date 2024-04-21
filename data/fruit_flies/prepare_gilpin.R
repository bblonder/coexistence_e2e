library(dplyr)
library(tidyr)

reformat_pairwise <- function(pairwise)
{
  # check symmetry of coexistence (*)
  test_coexistence = (pairwise=="*") == t(pairwise=="*")
  print(which(test_coexistence!=TRUE, arr.ind = TRUE))
  
  # check symmetry of competitive hierarchy
  test_hierarchy = (pairwise=="*") | ( (pairwise=="1") & t(pairwise=="0") ) | (pairwise=="0") & t(pairwise=="1")
  print(which(test_hierarchy!=TRUE, arr.ind=TRUE))
  
  test_extinct = (pairwise=="0") & t(pairwise=="0")
  print(which(test_extinct==TRUE, arr.ind=TRUE))
  
  # keep only upper triangular because we know it is symmetric
  pairwise[lower.tri(pairwise)] = NA
  
  pairwise_long = pairwise %>% 
    mutate(species_1 = row.names(pairwise)) %>%
    pivot_longer(-species_1,
                 names_to="species_2",
                 values_to="outcome") %>%
    filter(!is.na(outcome))
  
  return(pairwise_long)
}

pairwise_a = read.csv('figure_2.2a.csv',na=c("", " "), row.names=1, colClasses='character')
pairwise_a_long = reformat_pairwise(pairwise_a) %>%
  mutate(food='thick', temperature='25')

pairwise_b = read.csv('figure_2.2b.csv',na=c("", " "), row.names=1, colClasses='character')
pairwise_b_long = reformat_pairwise(pairwise_b) %>%
  mutate(food='thick', temperature='19')

pairwise_c = read.csv('figure_2.2c.csv',na=c("", " "), row.names=1, colClasses='character')
pairwise_c_long = reformat_pairwise(pairwise_c) %>%
  mutate(food='thin', temperature='25')

pairwise_d = read.csv('figure_2.2d.csv',na=c("", " "), row.names=1, colClasses='character')
pairwise_d_long = reformat_pairwise(pairwise_d) %>%
  mutate(food='thin', temperature='19')

pairwise_long_all = rbind(pairwise_a_long, pairwise_b_long, pairwise_c_long, pairwise_d_long)



# load in multispecies data
multispecies = read.csv('table_2.5.csv') %>%
  mutate(food='thick', temperature='25')



# make a combined data frame
species_all = names(pairwise_a) # pull from the most complete name list

df_empty = data.frame(matrix(data=0, nrow=1, ncol=length(species_all)+2+length(species_all)))
colnames(df_empty) = c(paste(species_all,"action",sep="."),'food.initial','temperature.initial',paste(species_all,"star",sep="."))

# make final data frame
result = NULL

# add in the pairwise experiments
for (i in 1:nrow(pairwise_long_all))
{
  sp1_this = pairwise_long_all[i,"species_1",drop = TRUE]
  sp2_this = pairwise_long_all[i,"species_2",drop = TRUE]
  outcome_this = pairwise_long_all[i,"outcome",drop = TRUE]
  food_this = pairwise_long_all[i,"food",drop = TRUE]
  temperature_this = pairwise_long_all[i,"temperature",drop = TRUE]
  
  df_this = df_empty
  df_this[1,paste(sp1_this,"action",sep=".")] = 0.5 # relative abundance
  df_this[1,paste(sp2_this,"action",sep=".")] = 0.5 # relative abundance 

  if (outcome_this == "*")
  {
    df_this[1,paste(sp1_this,"star",sep=".")] = 0.5 # scale this as a relative abundance, assume equal population sizes (no data available to predict more quantitatively)
    df_this[1,paste(sp2_this,"star",sep=".")] = 0.5 
  }
  else if (outcome_this == "1")
  {
    df_this[1,paste(sp1_this,"star",sep=".")] = 1
    df_this[1,paste(sp2_this,"star",sep=".")] = 0   
  }
  else if (outcome_this == "0")
  {
    df_this[1,paste(sp1_this,"star",sep=".")] = 0
    df_this[1,paste(sp2_this,"star",sep=".")] = 1   
  }
  else
  {
    stop('this combination should not exist')
  }
  
  df_this[1,"food.initial"] = food_this
  df_this[1,"temperature.initial"] = temperature_this
  
  
  result = rbind(result, df_this)
}

# make sure we have the competitive dominance 1s/0s correct - verified on 4/16/24 as row i excludes column j

# add in the multispecies as relative abundances
names_multispecies = multispecies %>% select(Merc:Gib) %>% names

result_multispecies = NULL

for (i in 1:nrow(multispecies))
{
  df_this = df_empty
  
  for (name_this in names_multispecies)
  {
    df_this[1, paste(name_this,"action",sep=".")] = multispecies[i,name_this] / 100 # scale as relative abundnace
    df_this[1, paste(name_this,"star",sep=".")] = multispecies[i,paste(name_this,"outcome",sep=".")] / 100 # this scales to relative abundance
  }
  
  df_this[1,"food.initial"] = multispecies[i,"food"]
  df_this[1,"temperature.initial"] = multispecies[i,"temperature"]
  
  result_multispecies = rbind(result_multispecies, df_this)
}


# put the names back in the expected letter format
result_all = rbind(result, result_multispecies)

write.csv(result_all, file='data_fruit_flies.csv', row.names=FALSE)
