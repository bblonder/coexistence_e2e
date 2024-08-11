
get_dataset_stats <- function(fn) 
{
  data = read.csv(fn)
  
  n_cases = nrow(data)
  
  n_cases_trimmed = quantile_max_trim(data) %>% 
    select(contains("outcome")) %>%
    na.omit %>%
    nrow
  
  n_na = n_cases - n_cases_trimmed
  
  n_sp = data %>% 
    select(contains("outcome")) %>% 
    ncol
  
  n_env = data %>% 
    select(contains("initial")) %>% 
    ncol
  
  n_env_levels = data %>% 
    select(contains("initial")) %>% 
    unique %>% 
    nrow
  
  n_combos = data %>%
    select(contains(c('action','initial'))) %>%
    unique %>% 
    nrow
  
  q_995 = data %>% 
    select(contains("outcome")) %>% 
    as.matrix %>%
    as.numeric %>%
    ifelse(.>1e7,NA,.) %>%
    quantile(0.995,na.rm=T) %>%
    as.numeric
  
  return(data.frame(fn, name = gsub('data_','',gsub('\\.csv','',basename(fn))),
                    n_cases, n_na, n_sp, n_env, n_env_levels, n_combos, q_995))
}
