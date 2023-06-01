freq_weight <- function(x)
{
  x_cut = cut(x, breaks=10)
  w = 1/table(x_cut)
  w = as.numeric(w)
  w[!is.finite(w)] = NA
  w <- w/sum(w,na.rm=T)
  names(w) = levels(x_cut)
  
  weights_by_x = as.numeric(w[x_cut])
  
  return(weights_by_x)
}