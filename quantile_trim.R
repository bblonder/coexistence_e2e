library(dplyr)

quantile_trim <- function(df, q=0.5*1e-2)
{
  quantiles = df %>%
    select(contains("star")) %>%
    as.matrix %>%
    as.numeric %>%
    quantile(c(q, 1 - q))
  
  print(quantiles)
  
  df_processed = df %>%
    mutate(across(contains("star"), function(x) { ifelse(x < quantiles[1] | x > quantiles[2], NA, x)}))
  
  num_nas_before = sum(is.na(df))
  print(num_nas_before)
  
  num_nas_processed = sum(is.na(df_processed))
  print(num_nas_processed)
  
  return(df_processed)
}
