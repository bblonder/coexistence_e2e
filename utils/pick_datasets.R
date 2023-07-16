pick_datasets <- function(name=NULL, rep=NULL, num_train=NULL, method=NULL, experimental_design=NULL, response_var=NULL)
{
  possible_files = dir(sprintf('outputs/statistical/%s',name),pattern='index*',full.names = TRUE)
  
  if(!is.null(rep))
  {
    ids = grep(pattern=sprintf('rep=%d_',rep),possible_files)
    possible_files = possible_files[ids]
  }
  
  if (!is.null(num_train))
  {
    ids = grep(pattern=sprintf('num_train=%d',num_train),possible_files)
    possible_files = possible_files[ids]
  }
  
  if (!is.null(method))
  {
    ids = grep(pattern=sprintf('method=%s_',method),possible_files)
    possible_files = possible_files[ids]
  }

  if (!is.null(experimental_design))
  {
    ids = grep(pattern=sprintf('experimental_design=%s',experimental_design),possible_files)
    possible_files = possible_files[ids]
  }
  
  if (!is.null(response_var))
  {
    ids = grep(pattern=sprintf('response=%s',response_var),possible_files)
    possible_files = possible_files[ids]
  }
  
  if (length(possible_files)==0)
  {
    return(NULL)
  }
  else
  {
    files_obs_train = possible_files[grep("output=obs_train",possible_files)]
    files_pred_train = possible_files[grep("output=pred_train",possible_files)]
    files_exp_train = possible_files[grep("experiment_train",possible_files)]
    
    files_obs_test = possible_files[grep("output=obs_test",possible_files)]
    files_pred_test = possible_files[grep("output=pred_test",possible_files)]
    files_exp_test = possible_files[grep("experiment_test",possible_files)]
    
    return(data.frame(files_obs_train = files_obs_train,
                      files_pred_train = files_pred_train,
                      files_exp_train = files_exp_train,
                      files_obs_test = files_obs_test,
                      files_pred_test = files_pred_test,
                      files_exp_test = files_exp_test,
                      name = name))
  }
}
