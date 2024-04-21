library(reshape)
library(dplyr)
library(readxl)

# data are from Friedman and Gore

sheets_mono = excel_sheets("data/monoculture_timeSeries.xlsx")
data_singlet = do.call("rbind",lapply(1:length(sheets_mono), function(i) {
  abund_final = read_excel('data/monoculture_timeSeries.xlsx',sheet=i) %>% 
    select(A:H) %>%
    tail(1) %>%
    t
  
  xoutcome = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  xoutcome[,i] = abund_final
  xoutcome = data.frame(xoutcome)
  names(xoutcome) = paste(sheets_mono,"outcome",sep=".")
  
  x = matrix(0, nrow=nrow(xoutcome),ncol=length(sheets_mono))
  x[,i] = 1
  x = data.frame(x)
  names(x) = sheets_mono
  
  return(data.frame(x, 
                    xoutcome))
  }))


sheets_pair = excel_sheets('data/pair_timeSeries.xlsx')

data_pair = do.call("rbind",lapply(1:length(sheets_pair), function(i) {
  abund_final = read_excel('data/pair_timeSeries.xlsx',sheet=i)
  names(abund_final)[1] = 'time'
  abund_final = abund_final %>% 
    filter(time==5) %>% # get ifnal times
    select(-time)
  
  xoutcome = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  xoutcome = data.frame(xoutcome)
  names(xoutcome) = sheets_mono
  xoutcome[,names(abund_final)] = abund_final
  names(xoutcome) = paste(names(xoutcome),"outcome",sep=".")
  
  x = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  x = data.frame(x)
  names(x) = sheets_mono
  
  x[,names(abund_final)] = 1
  
  return(data.frame(x, 
                    xoutcome))
  }))



sheets_triplet = excel_sheets('data/trio_lastTransfer.xlsx')

data_triplet = do.call("rbind",lapply(1:length(sheets_triplet), function(i) {
  abund_final = read_excel('data/trio_lastTransfer.xlsx',sheet=i)
  names(abund_final)[1] = 'time'
  abund_final = abund_final %>% 
    filter(time==5) %>% # get ifnal times
    select(-time)
  
  xoutcome = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  xoutcome = data.frame(xoutcome)
  names(xoutcome) = sheets_mono
  xoutcome[,names(abund_final)] = abund_final
  names(xoutcome) = paste(names(xoutcome),"outcome",sep=".")
  
  x = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  x = data.frame(x)
  names(x) = sheets_mono
  
  x[,names(abund_final)] = 1
  
  return(data.frame(x, 
                    xoutcome))
}))


sheets_all_and_dropouts = excel_sheets('data/7and8Species_lastTransfer.xlsx')
data_all_and_dropouts = do.call("rbind",lapply(1:length(sheets_all_and_dropouts), function(i) {
  abund_final = read_excel('data/7and8Species_lastTransfer.xlsx',sheet=i)
  names(abund_final)[1] = 'time'
  abund_final = abund_final %>% 
    filter(time==5) %>% # get ifnal times
    select(-time)
  
  xoutcome = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  xoutcome = data.frame(xoutcome)
  names(xoutcome) = sheets_mono
  xoutcome[,names(abund_final)] = abund_final
  names(xoutcome) = paste(names(xoutcome),"outcome",sep=".")
  
  x = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  x = data.frame(x)
  names(x) = sheets_mono
  
  x[,names(abund_final)] = 1
  
  return(data.frame(x, 
                    xoutcome))
}))




data_combined = rbind(data_singlet, data_pair, data_triplet, data_all_and_dropouts)

for (i in grep("outcome",names(data_combined)))
{
  data_combined[,i] = as.numeric(data_combined[,i])
}

# remove empty cases
which_rows_na = data_combined %>%
  select(contains("outcome")) %>% 
  rowSums %>%
  is.na %>%
  which
data_combined = data_combined[-which_rows_na,]

# rename action columns
num_species = length(grep("outcome$",names(data_combined)))
names(data_combined)[1:num_species] = paste(names(data_combined)[1:num_species],"action",sep=".")

# add environment dummy column
data_combined = data_combined %>%
  mutate(environment.initial = 0) %>%
  select(contains("action"),contains("initial"),contains("outcome"))

write.csv(data_combined, file='data_soil_bacteria.csv',row.names=F)
