# try to visualize interactions and explanations
source('src/coexistence_love.R')
source('src/configs.R')

num_species_M_11 = 11
num_train_M_11 = 500
data_assemblages_M_11 = read.csv('data/glv/assemblages_M_11.csv')
data_assemblages_M_11 = clean_input_data(data_assemblages_M_11)
data_assemblages_M_11 = get_state_assemblages_mapping(num_species, data_assemblages_M_11)

idx_train = generate_state_idxs_train(experimental_design = 'mixed',
                                      num_train = num_train_M_11,
                                      assemblages = data_assemblages_M_11,
                                      method = 'rf',
                                      num_species = num_species_M_11,
                                      hyperparams = MODEL_HYPERPARAMS)

z = perform_prediction_experiment_single(
  predictor_variable = '_abundance',
  assemblages = data_assemblages_M_11,
  num_train = num_train_M_11,
  state_idxs_train = idx_train,
  state_idxs_test = idx_train,
  method = 'rf',
  num_species = num_species_M_11)

w = find.interaction(z$model, 
                     trace=TRUE,
                     verbose=TRUE,
                     sorted=FALSE,
                     method='vimp')

w_mat = matrix(NA,
               nrow=num_species_M_11,
               ncol=num_species_M_11)
w_mat[upper.tri(w_mat,diag = FALSE)] = w[,"Difference"] 



r_M = c(0.36807,0.31023,0.3561,0.54006,0.70898,0.47064,0.2297,0.83005,0.39181,0.29075,0.32367)
names(r_M) = letters[1:length(r_M)]
A_M = c(-0.20516,0.098398,0.16739,-0.16461,-0.14341,0.019881,-0.51535,-0.39162,0.34635,0.0088853,-0.26894,0.062123,-0.10489,-0.043011,-0.15466,-0.1872,0.027031,-0.45919,-0.41388,0.3013,0.022081,-0.19657,0.14373,-0.19203,-0.10162,-0.13971,-0.16537,0.013651,-0.50414,-0.7724,0.29257,-0.005959,-0.20645,0.22403,0.13813,0.00045883,-0.83125,-0.2238,0.22027,-0.20529,-1.0097,0.66639,-0.038986,-0.40032,-0.18016,-0.051261,-5.03E-05,-0.054212,-0.70858,0.016198,-0.50756,0.55363,0.15757,0.22438,0.10635,-0.11159,-0.03721,-0.042591,0.041044,0.26134,-0.42266,-0.18536,-0.43231,0.1647,-0.061038,-0.26461,-0.12669,-0.18576,-0.12222,0.3809,0.4003,-0.16078,-1.2124,1.3897,-0.37922,0.19189,-0.096352,-0.071257,0.00060448,0.080355,-0.4548,-0.50349,0.16899,-0.56222,-4.3508,0.44315,-0.22341,-0.2074,-0.037541,-0.033333,-0.049912,-0.090424,-0.10211,0.03229,-0.18179,-0.30301,-0.055765,0.01436,-0.0076697,-0.04225,-0.013105,0.02398,-0.11784,-0.32893,0.020748,0.054767,-2.0963,0.11124,-0.19213,0.023816,   -0.3742,0.27843,0.24887,-0.16829,0.08399,0.033691,-0.23242,-0.39513,0.31454,-0.038764,-0.3841)
A_M = matrix(A_M, nrow=length(r_M),ncol=length(r_M),dimnames=list(letters[1:length(r_M)],letters[1:length(r_M)]))

plot(w_mat,abs(A_M))
plot(w_mat,abs(t(A_M)))


# q = partial.rfsrc(z$model,
#               partial.xvar='a')
