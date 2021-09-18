library(dplyr)
library(tidyr)
library(progress)
library(parallel)
library(data.table)

# we assume the row and column orders are the same
params_vital = read.csv('speciesvitalrates.csv') %>% 
					select(-species)
params_alpha = read.csv('alpha_estimates_row_is_target.csv') %>%
					select(-X) %>%
					as.matrix
params_alpha_gapfill = params_alpha
params_alpha_gapfill[is.na(params_alpha_gapfill)] = mean(params_alpha,na.rm=T)

params_godoy = list(vital=params_vital, alpha=params_alpha_gapfill)


# make assemblages
generate_assemblages <- function(n, labels=letters)
{
  df = expand.grid(replicate(n, 0:1, simplify = FALSE))
  names(df) = labels[1:n]
  
  x_star = as.data.frame(matrix(data=0,nrow=nrow(df),ncol=ncol(df)))
  names(x_star) = paste(names(df),"star",sep=".")
  
  df_final = data.frame(df, 
                        stable=NA, 
                        feasible=NA,
                        richness=NA,
                        x_star)
    
  return(df_final)
}


assign_params <- function(assemblage, params)
{
  stopifnot(nrow(assemblage)==1)
  # find species that are present
  species_indices_present = which(as.numeric(assemblage)==1)

  # pick subset of parameters (assuming that the parameters don't change when subsetting)
  vital_this = params$vital[species_indices_present, , drop=FALSE]
  alpha_this = params$alpha[species_indices_present, species_indices_present, drop=FALSE]
  
  return(list(vital= vital_this,alpha= alpha_this))
}

do_simulation_annual_plant <- function(vital, alpha, nstep=200)
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


fill_in_assemblages <- function(assemblages, params)
{
  #pb <- progress_bar$new(format='[:bar] :elapsed :eta :percent :spin')
  
  n = log2(nrow(assemblages))
  result = mclapply(1:nrow(assemblages), function(i)
  {
    cat('.\n')
    params_this_row = assign_params(assemblage = assemblages[i,1:n], params = params)
    
    names_this = letters[1:n]
    names_this = names_this[which(as.numeric(assemblages[i,1:n])==1)]
    
    nt = do_simulation_annual_plant(vital = params_this_row$vital, alpha = params_this_row$alpha)
    dimnames(nt) = list(NULL, names_this)
    
    x_star = nt[nrow(nt),,drop=TRUE]

    assemblages[i,"stable"] = mean(apply(tail(nt,5), 2, function(x) {sd(x) / mean(x)}  )) < 0.5
    assemblages[i,"feasible"] = all(nt[nrow(nt),] >= 0)
    
    if (length(x_star) > 0)
    {
      names(x_star) = paste(names(x_star),"star",sep=".")
      
      assemblages[i,names(x_star)] = x_star
    }
    
    #pb$update(i/nrow(assemblages))
    
    assemblages[i,"richness"] = sum(x_star > 0.01)
    
    return(assemblages[i,])
  }, mc.cores=16)
  result_final = rbindlist(result)
  #pb$terminate()
  
  return(result_final)
}

assemblages_godoy = generate_assemblages(n = nrow(params_vital))

assemblages_godoy = fill_in_assemblages(params = params_godoy, assemblages = assemblages_godoy)

write.csv(assemblages_godoy,file='assemblages_annual_plant_18.csv',row.names=FALSE)

