library(dplyr)
library(tidyr)
library(deSolve)
library(ggplot2)
library(reshape2)

generate_assemblages <- function(n, labels=letters)
{
  df = expand.grid(replicate(n, 0:1, simplify = FALSE))
  names(df) = labels[1:n]
  
  x_outcome = as.data.frame(matrix(data=0,nrow=nrow(df),ncol=ncol(df)))
  names(x_outcome) = paste(names(df),"outcome",sep=".")
  
  df_final = data.frame(df, 
                        stable=NA, 
                        feasible=NA,
                        richness=NA,
                        x_outcome)
    
  return(df_final)
}

generate_params <- function(n, labels=letters, A.diag=-1,
                            distribution.A='norm',
                            A.norm.mean=-0.2, A.norm.sd=0.2,  
                            A.lnorm_neg.meanlog=log(0.25), A.lnorm_neg.sdlog=log(1.5),
                            A.lnorm_pos.meanlog=log(0.25), A.lnorm_pos.sdlog=log(1.5),
                            distribution.r='norm',
                            r.norm.mean=1, r.norm.sd=0.5, 
                            r.lnorm_pos.meanlog=log(1), 
                            r.lnorm_pos.sdlog=log(2))
{
  # interaction matrix
  A = matrix(data=NA,nrow=n,ncol=n)
  dimnames(A) = list(labels[1:n],labels[1:n])
  
  # intrinsic growth rate vector
  r = rep(NA, n)
  names(r) = labels[1:n]
  
  # fill A
  if (distribution.A=='norm')
  {
    A[] = rnorm(n=n*n,mean=A.norm.mean,sd=A.norm.sd) 
  }
  else if (distribution.A=='lnorm_neg')
  {
    A[] = -1*rlnorm(n=n*n,meanlog=A.lnorm_neg.meanlog,sdlog=A.lnorm_neg.sdlog) 
  }
  else if (distribution.A=='lnorm_pos')
  {
    A[] = rlnorm(n=n*n,meanlog=A.lnorm_pos.meanlog,sdlog=A.lnorm_pos.sdlog)
  }
  else
  {
    stop('distribution not found')
  }
  
  # fill r
  if (distribution.r=='norm')
  {
    r[] = rnorm(n=n,mean=r.norm.mean,sd=r.norm.sd)
  }
  else if (distribution.r=='lnorm_pos')
  {
    r[] = rlnorm(n=n,meanlog=r.lnorm_pos.meanlog,sdlog=r.lnorm_pos.sdlog)
  }
  else
  {
    stop('distribution not found')
  }
  
  # set diagonal for A if requested
  if (!is.null(A.diag))
  {
    diag(A) <- rep(A.diag,n)
  }
 
  # return combo
  return(list(A=A,r=r))
}


assign_params <- function(assemblage, params)
{
  stopifnot(nrow(assemblage)==1)
  # find species that are present
  species_indices_present = which(as.numeric(assemblage)==1)
  
  # pick subset of parameters (assuming that the parameters don't change when subsetting)
  A_this = params$A[species_indices_present, species_indices_present, drop=FALSE]
  r_this = params$r[species_indices_present]
  
  return(list(A=A_this,r=r_this))
}

determine_fixed_point <- function(params)
{
  if (length(params$r) > 0)
  {
    # do forward simulation
    LVmod = function(t, N, params) {
        r = params$r
        A = params$A
        dNdt = N * (r + A %*% N)
        list(dNdt)
    } 
    
    time = seq(from=0, to=1000, by=0.1)
    x_0 = rep(1, length(params$r))
    
    simulation = ode(
        y = x_0,
        times = time,
        func = LVmod,
        parms = params
      ) %>%
      as.data.frame %>%
      as_tibble
    
    # sim_to_plot = melt(simulation,id.vars='time')
    # print(str(sim_to_plot))
    # g = ggplot(sim_to_plot,aes(x=time,y=value,color=variable)) +
    #   geom_line()
    # ggsave(g,file=sprintf('~/Downloads/test_%f.png',runif(1)))
    
    x_outcome =
      simulation %>%
      filter(time == max(time)) %>%
      select(-time) %>%
      as.numeric()
    
    names(x_outcome) = paste(names(params$r),"outcome",sep=".")
  }
  else
  {
      x_outcome = NULL
  }

  return(x_outcome)
}

determine_feasibility <- function(params, x_outcome)
{
  if (length(params$r) > 0)
  {
    feasibility = all(x_outcome > 0)
  }
  else
  {
    feasibility = TRUE
  }
  
  return(feasibility)
}

determine_stability <- function(params, x_outcome)
{
  if (length(params$r) > 0)
  {
    Jacobian = diag(x=x_outcome,nrow=nrow(params$A),ncol=ncol(params$A)) * params$A
    lambda = eigen(Jacobian)$values
    
    stability = all(Re(lambda) < 0)
  }
  else
  {
    stability = TRUE
  }
  
  return(stability)
}

fill_in_assemblages <- function(assemblages, params)
{
  n = log2(nrow(assemblages))
  for (i in 1:nrow(assemblages)) # skip the no-species scenario
  {
    message(i/nrow(assemblages))
    params_this_row = assign_params(assemblage = assemblages[i,1:n], params = params)
    
    x_outcome = determine_fixed_point(params_this_row)
    
    assemblages[i,"stable"] = determine_stability(params_this_row, x_outcome=x_outcome)
    assemblages[i,"feasible"] = determine_feasibility(params_this_row, x_outcome=x_outcome)
    
    if (!is.null(x_outcome))
    {
      assemblages[i,names(x_outcome)] = x_outcome
    }
    
    # set abundance of 0.01 as the minimum threshold to 'count' for richness
    assemblages[i,"richness"] = sum(x_outcome > 0.01)
  }
  return(assemblages)
}


# 
# 
# set.seed(1) # replicability
# nsp_glv_16 <- 16
# warning('change # of species back!')
# params_glv_16 = generate_params(n = nsp_glv_16, A.norm.mean = -0.2, r.norm.mean = 1.5)
# assemblages_glv_16 = fill_in_assemblages(params = params_glv_16, assemblages = generate_assemblages(n = nsp_glv_16))
# write.csv(assemblages_glv_16,file='assemblages_glv_16.csv',row.names=FALSE)
# # check counts
# assemblages_glv_16 %>% select(richness, stable) %>% table
# 
# 


process_for_love <- function(df, names)
{
  df = df %>% select(-stable, -feasible, -richness)
  stopifnot(ncol(df)==2*length(names))
  names(df)[1:length(names)] = paste(names, "action",sep=".")
  names(df)[(length(names)+1):ncol(df)] = paste(names, "outcome",sep=".")
  
  df = df %>%
    mutate(environment.initial = 0) %>%
    select(contains("action"),contains("initial"),contains("outcome"))
  
  return(df)
}



# 12 Member Human Gut Community from Venturelli et al. 
# https://www.embopress.org/doi/full/10.15252/msb.20178157
r_H = c(0.2453, 0.4777, 0.5979, 0.5841, 0.4573, 0.2464, 0.5025, 0.2321, 0.4021, 0.1558, 0.2192, 0.2375)
names(r_H) = letters[1:length(r_H)]
A_H = c(-0.9118,-0.2145,-0.2718,-0.2275,-0.1294,-0.3058,-0.3478,-0.9002,0.1764,-0.5455,-0.2307,-0.5286,0,-0.7339,-0.6235,-0.6317,-0.515,0,-0.5069,0,1.7566,0,-0.2087,0,0,-0.8193,-0.9067,-0.7538,-0.7552,0,-0.0864,0,2.2719,-0.7375,-0.7033,0,0,-0.9208,-0.8155,-0.8804,-0.5498,0,0.0656,0,3.3782,0,-0.7822,0,0.1366,-0.5556,-0.6416,-0.5837,-0.6597,-0.6566,-0.0468,-0.1062,1.3026,0,-0.6379,0,0.4526,-0.2776,-0.2736,-0.2614,-0.1679,-0.829,0.3118,0,-0.4475,0,-1.121,-0.6709,0,-0.4645,-0.6319,-0.1511,0.039,-0.2413,-1.4543,-2.157,0.0239,-0.4394,-0.5069,-0.771,0,-0.2028,-0.1999,-0.1763,-0.026,0.042,-0.151,-1.2535,0.1756,0,0,-0.4333,0.6924,-0.1141,-0.1688,-0.1241,-0.0493,-1.098,0,-0.4084,-2.4418,0,-0.1531,-1.0774,1.3433,-0.0203,0,-0.0613,-0.0296,0,1.0831,0,-0.1385,-1.2705,-0.1676,0,0.9613,-0.0993,-0.0685,0.2313,0.7585,0,0.4481,1.0147,-0.768,0,-1.0382,-0.402,0,-0.2646,-0.3033,-0.3237,-0.2017,-0.5596,0.265,-0.9771,-0.9041,-0.8171,-0.4053,-0.6217)
A_H = matrix(A_H, nrow=length(r_H),ncol=length(r_H),dimnames=list(letters[1:length(r_H)],letters[1:length(r_H)]))
params_H = list(A=A_H,r=r_H)
assemblages_H = fill_in_assemblages(params = params_H, assemblages = generate_assemblages(n = length(r_H)))

names_H = c("BH","CA","BU","PC","BO","BV","BT","EL","FP","CH","DP","ER")

assemblages_H_final = process_for_love(assemblages_H, names_H)

write.csv(assemblages_H_final,file='data_human_gut.csv',row.names=FALSE)


# 11 Member Mouse Gut Community from MDSINE series of papers, Bucci et al. 
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003388
r_M = c(0.36807,0.31023,0.3561,0.54006,0.70898,0.47064,0.2297,0.83005,0.39181,0.29075,0.32367)
names(r_M) = letters[1:length(r_M)]
A_M = c(-0.20516,0.098398,0.16739,-0.16461,-0.14341,0.019881,-0.51535,-0.39162,0.34635,0.0088853,-0.26894,0.062123,-0.10489,-0.043011,-0.15466,-0.1872,0.027031,-0.45919,-0.41388,0.3013,0.022081,-0.19657,0.14373,-0.19203,-0.10162,-0.13971,-0.16537,0.013651,-0.50414,-0.7724,0.29257,-0.005959,-0.20645,0.22403,0.13813,0.00045883,-0.83125,-0.2238,0.22027,-0.20529,-1.0097,0.66639,-0.038986,-0.40032,-0.18016,-0.051261,-5.03E-05,-0.054212,-0.70858,0.016198,-0.50756,0.55363,0.15757,0.22438,0.10635,-0.11159,-0.03721,-0.042591,0.041044,0.26134,-0.42266,-0.18536,-0.43231,0.1647,-0.061038,-0.26461,-0.12669,-0.18576,-0.12222,0.3809,0.4003,-0.16078,-1.2124,1.3897,-0.37922,0.19189,-0.096352,-0.071257,0.00060448,0.080355,-0.4548,-0.50349,0.16899,-0.56222,-4.3508,0.44315,-0.22341,-0.2074,-0.037541,-0.033333,-0.049912,-0.090424,-0.10211,0.03229,-0.18179,-0.30301,-0.055765,0.01436,-0.0076697,-0.04225,-0.013105,0.02398,-0.11784,-0.32893,0.020748,0.054767,-2.0963,0.11124,-0.19213,0.023816,   -0.3742,0.27843,0.24887,-0.16829,0.08399,0.033691,-0.23242,-0.39513,0.31454,-0.038764,-0.3841)
A_M = matrix(A_M, nrow=length(r_M),ncol=length(r_M),dimnames=list(letters[1:length(r_M)],letters[1:length(r_M)]))
params_M = list(A=A_M,r=r_M)
assemblages_M = fill_in_assemblages(params = params_M, assemblages = generate_assemblages(n = length(r_M)))

names_M = c("Bar","undLac","uncLac","Oth","Bla","undMol","Akk","Cop","Clodif","Ent","undEnt")

assemblages_M_final = process_for_love(assemblages_M, names_M)


write.csv(assemblages_M_final,file='data_mouse_gut.csv',row.names=FALSE)

