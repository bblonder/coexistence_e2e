library(dplyr)

quantile_max_trim <- function(df, q=0.5*1e-2, max_val=1e7, round_zero = TRUE)
{
  cols_star = grep("star",names(df))
  for (i in cols_star)
  {
    df[df[,i] < 0 & df[,i] > -1e-6, i] = 0
    df[df[,i] > 1e7, i] = NA
  }
  
  quantiles = df %>%
    select(contains("star")) %>%
    as.matrix %>%
    as.numeric %>%
    quantile(c(q, 1 - q), na.rm=T)
  
  df_processed = df %>%
    mutate(across(contains("star"), function(x) { ifelse(x < quantiles[1] | x > quantiles[2], NA, x)}))
  
  num_nas_before = sum(is.na(df))
  #print(num_nas_before)
  
  num_nas_processed = sum(is.na(df_processed))
  #print(num_nas_processed)
  
  return(df_processed)
}
