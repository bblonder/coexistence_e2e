num_species = 12
assemblages = read.csv('data/glv/assemblages_H_12.csv')
training_rows = 1:nrow(assemblages)
# training_rows = get_single_species_and_leave_one_out_rows(assemblages, num_species)

initial_cols = 1:num_species
existence_cols = names(assemblages)[grep("star", names(assemblages))]

initial_vec = as.numeric(assemblages_row[initial_cols])
glv_training = create_glv_dataset(
  assemblages[training_rows,], num_species)
X = data.matrix(glv_training)
y = data.matrix(rep(-1, nrow(glv_training)))
a_eff = as.numeric(ginv(X) %*% y)
a_eff_mat = matrix(a_eff, nrow = num_species, byrow=TRUE)

test_row = assemblages[420,]
test_initial = as.numeric(test_row[,initial_cols])
test_existence = as.numeric(test_row[,existence_cols])
initial_idx = which(test_initial == 1)
a_sub_mat = a_eff_mat[initial_idx,initial_idx]

test_existence_pred_raw = as.numeric(-ginv(a_sub_mat) %*% rep(1, length(initial_idx)))
test_existence_pred = rep(0, num_species)
test_existence_pred[initial_idx] = clamp(test_existence_pred_raw, lower = 0)
test_existence
test_existence_pred

####
num_species = 8
input_file = read.csv('data/friedman_gore/data_friedman_gore.csv')
assemblages = clean_input_data(input_file)
assemblages = get_state_assemblages_mapping(num_species, assemblages)
assemblages = assemblages %>% 
  slice_sample(n = 2, by = state_idx)
training_state_idxs = 1:2^num_species
predict_state_idxs = 1:2^num_species
assemblages_training = data.frame(
  get_assemblages_subset_from_state_idxs(
    training_state_idxs, assemblages))

glv_wrapper = fit_glv_baseline(assemblages, training_state_idxs, num_species)
glv_predictions = predict_glv("_abundance", glv_wrapper, assemblages, predict_state_idxs, num_species)
assemblages_training[diff_columns] = assemblages_training[star_columns] - glv_predictions


#####