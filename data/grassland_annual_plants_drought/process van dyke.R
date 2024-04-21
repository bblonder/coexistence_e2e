library(dplyr)
library(tidyr)
library(progress)
library(data.table)
library(ggplot2)

# data are from https://zenodo.org/records/7083314

data_s_g_all = read.csv('s_g_data.csv', row.names = 1)

data_alpha_all = read.csv('nls_boot_pairs_1000_full_model.csv') %>%
  select(focal, competitor, treatment, alpha, lambda)

prep_parameters <- function(data_s_g, data_alpha, treatment_this)
{
  alpha_matrix_treatment_this = data_alpha %>% 
    filter(treatment==treatment_this) %>% 
    select(-treatment, -lambda) %>% 
    pivot_wider(names_from=focal, values_from=alpha) %>%
    as.data.frame
  row.names(alpha_matrix_treatment_this) = alpha_matrix_treatment_this$competitor
  alpha_matrix_treatment_this = alpha_matrix_treatment_this %>%
    select(-competitor) %>%
    as.matrix
  
  lambda_this = data_alpha %>%
    select(focal, treatment, lambda) %>%
    filter(treatment==treatment_this) %>%
    unique %>%
    select(-treatment)
  vital_this = data_s_g %>% 
    left_join(lambda_this,by='focal')
  row.names(vital_this) = vital_this$focal
  vital_this = vital_this %>%
    select(-focal)
  
  return(list(vital=vital_this, alpha=alpha_matrix_treatment_this, treatment=treatment_this))
}



params_treatment_1 = prep_parameters(data_s_g_all, data_alpha_all, treatment=1)
params_treatment_2 = prep_parameters(data_s_g_all, data_alpha_all, treatment=2)




names_species = data_s_g_all$focal





# make assemblages
generate_assemblages <- function(labels)
{
  n = length(labels)
  
  df = expand.grid(replicate(n, 0:1, simplify = FALSE))
  names(df) = labels
  
  df_outcome = as.data.frame(matrix(data=0,nrow=nrow(df),ncol=ncol(df)))
  names(df_outcome) = paste(names(df),"outcome",sep=".")
  
  df_final = data.frame(df, 
                        df_outcome)
  
  return(df_final)
}


assign_params <- function(assemblage, params)
{
  stopifnot(nrow(assemblage)==1)
  # find species that are present
  species_indices_present = which(as.numeric(assemblage)==1)
  names_species_present = names(assemblage)[species_indices_present]
  print(names_species_present)
  
  # pick subset of parameters (assuming that the parameters don't change when subsetting)
  vital_this = params$vital[species_indices_present, , drop=FALSE]
  alpha_this = params$alpha[species_indices_present, species_indices_present, drop=FALSE]
  
  return(list(vital= vital_this,alpha= alpha_this))
}






do_simulation_annual_plant <- function(vital, alpha, nstep=1000)
{
  nsp = nrow(vital)
  
  n = matrix(NA, nrow=nstep,ncol=nsp)
  
  if (nrow(vital) > 0)
  {
    # set initial abundance
    n[1,] = 1
    
    for (t in 1:(nrow(n)-1))
    {
      fec = rep(NA, nsp)
      for (i in 1:nsp)
      {
        term2 = 0
        for (j in 1:nsp)
        {
          term2 = term2 + alpha[i,j]*vital$g[j]*n[t,j]
        }
        fec[i] = vital$lambda[i] / ( 1 + term2  )
        n[t+1,i] = n[t,i] * ( (1 - vital$g[i])*vital$s[i] + vital$g[i]*fec[i] )
      }
    }
  }
  
  return(n)
}

plot_time_series <- function(nt, name)
{
  if (ncol(nt) > 0)
  {
    nt_long = nt %>% 
      as.data.frame %>%
      mutate(time=row_number()) %>%
      pivot_longer(-time)
    
    g = ggplot(nt_long, aes(x=time,y=value,color=name)) + 
      geom_line() +
      theme_bw() +
      ylim(0,8000)
    ggsave(g, file=sprintf('~/Downloads/%s.pdf',name))
  }
}


fill_in_assemblages <- function(assemblages, params, treatment_label)
{
  n = log2(nrow(assemblages))

  result = lapply(1:nrow(assemblages), function(i)
  {
    cat(sprintf('%d / %d\n',i, nrow(assemblages)))
    params_this_row = assign_params(assemblage = assemblages[i,1:n], params = params)
    
    names_this = names_species[which(as.numeric(assemblages[i,1:n])==1)]
    
    nt = do_simulation_annual_plant(vital = params_this_row$vital, alpha = params_this_row$alpha)
    dimnames(nt) = list(NULL, names_this)
    
    plot_time_series(nt, paste(c(names_this,treatment_label),collapse="."))
    
    x_outcome = nt[nrow(nt),,drop=TRUE]

    if (length(x_outcome) > 0)
    {
      names(x_outcome) = paste(names(x_outcome),"outcome",sep=".")
      
      assemblages[i,names(x_outcome)] = x_outcome
    }
    
    return(assemblages[i,])
  })
  result_final = rbindlist(result)
  
  names(result_final)[1:n] = paste(names_species,"action",sep=".")
  result_final = result_final %>%
    mutate(treatment.initial = treatment_label) %>%
    select(contains("action"),contains("initial"),contains("outcome"))

  return(result_final)
}

assemblages_treatment_1 = fill_in_assemblages(params = params_treatment_1, 
                                              assemblages = generate_assemblages(names_species),
                                              treatment_label='wet') # according to nls_orig_data.R 1==wet,2==dry
assemblages_treatment_2 = fill_in_assemblages(params = params_treatment_2, 
                                              assemblages = generate_assemblages(names_species),
                                              treatment_label='dry')

assemblages_all = rbind(assemblages_treatment_1, assemblages_treatment_2)

write.csv(assemblages_all,file='data_grassland_annual_plants_drought.csv',row.names=FALSE)

