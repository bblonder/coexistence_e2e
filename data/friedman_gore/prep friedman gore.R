library(reshape)
library(dplyr)
library(readxl)

sheets_mono = excel_sheets("data/monoculture_timeSeries.xlsx")
data_singlet = do.call("rbind",lapply(1:length(sheets_mono), function(i) {
  abund_final = read_excel('data/monoculture_timeSeries.xlsx',sheet=i) %>% 
    select(A:H) %>%
    tail(1) %>%
    t
  
  xstar = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  xstar[,i] = abund_final
  xstar = data.frame(xstar)
  names(xstar) = paste(sheets_mono,"star",sep=".")
  
  x = matrix(0, nrow=nrow(xstar),ncol=length(sheets_mono))
  x[,i] = 1
  x = data.frame(x)
  names(x) = sheets_mono
  
  return(data.frame(x, 
                    stable=NA, 
                    feasible=NA, 
                    richness=NA, 
                    xstar))
  }))


sheets_pair = excel_sheets('data/pair_timeSeries.xlsx')

data_pair = do.call("rbind",lapply(1:length(sheets_pair), function(i) {
  abund_final = read_excel('data/pair_timeSeries.xlsx',sheet=i)
  names(abund_final)[1] = 'time'
  abund_final = abund_final %>% 
    filter(time==5) %>% # get ifnal times
    select(-time)
  
  xstar = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  xstar = data.frame(xstar)
  names(xstar) = sheets_mono
  xstar[,names(abund_final)] = abund_final
  names(xstar) = paste(names(xstar),"star",sep=".")
  
  x = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  x = data.frame(x)
  names(x) = sheets_mono
  
  x[,names(abund_final)] = 1
  
  return(data.frame(x, 
                    stable=NA, 
                    feasible=NA, 
                    richness=NA, 
                    xstar))
  }))



sheets_triplet = excel_sheets('data/trio_lastTransfer.xlsx')

data_triplet = do.call("rbind",lapply(1:length(sheets_triplet), function(i) {
  abund_final = read_excel('data/trio_lastTransfer.xlsx',sheet=i)
  names(abund_final)[1] = 'time'
  abund_final = abund_final %>% 
    filter(time==5) %>% # get ifnal times
    select(-time)
  
  xstar = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  xstar = data.frame(xstar)
  names(xstar) = sheets_mono
  xstar[,names(abund_final)] = abund_final
  names(xstar) = paste(names(xstar),"star",sep=".")
  
  x = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  x = data.frame(x)
  names(x) = sheets_mono
  
  x[,names(abund_final)] = 1
  
  return(data.frame(x, 
                    stable=NA, 
                    feasible=NA, 
                    richness=NA, 
                    xstar))
}))


sheets_all_and_dropouts = excel_sheets('data/7and8Species_lastTransfer.xlsx')
data_all_and_dropouts = do.call("rbind",lapply(1:length(sheets_all_and_dropouts), function(i) {
  abund_final = read_excel('data/7and8Species_lastTransfer.xlsx',sheet=i)
  names(abund_final)[1] = 'time'
  abund_final = abund_final %>% 
    filter(time==5) %>% # get ifnal times
    select(-time)
  
  xstar = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  xstar = data.frame(xstar)
  names(xstar) = sheets_mono
  xstar[,names(abund_final)] = abund_final
  names(xstar) = paste(names(xstar),"star",sep=".")
  
  x = matrix(0, nrow=nrow(abund_final),ncol=length(sheets_mono))
  x = data.frame(x)
  names(x) = sheets_mono
  
  x[,names(abund_final)] = 1
  
  return(data.frame(x, 
                    stable=NA, 
                    feasible=NA, 
                    richness=NA, 
                    xstar))
}))




data_combined = rbind(data_singlet, data_pair, data_triplet, data_all_and_dropouts)
names(data_combined)[1:length(sheets_mono)] = letters[1:length(sheets_mono)]
names(data_combined)[grep("star",names(data_combined))] = paste(letters[1:length(sheets_mono)], "star",sep=".")

for (i in grep("star",names(data_combined)))
{
  data_combined[,i] = as.numeric(data_combined[,i])
}

data_combined = data_combined %>%
  mutate(richness = apply(data_combined[, grep("star",names(data_combined))], 1, function(x) { sum(x>0, na.rm=T)}))

# remove empty cases
which_rows_na = data_combined %>%
  select(contains("star")) %>% 
  rowSums %>%
  is.na %>%
  which
data_combined = data_combined[-which_rows_na,]

write.csv(data_combined, file='data_friedman_gore.csv',row.names=F)
