library(MASS)
library(raster)
library(ranger)
library(dplyr)
library(ggplot2)
library(tidyr)
library(caret)
library(randomForestSRC)
library(parallel)
library(data.table)
library(e1071)
library(vegan)
library(datastructures)
library(RANN)
library(conflicted)
conflict_prefer("union", "base")
conflict_prefer("intersect", "base")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Source utility functions
source('utils/freq_weight.R')
source('utils/skill_statistics.R')
source('utils/log_seq.R')
source('utils/quantile_trim.R')


get_predictor_columns <- function(
  predictor_variable,
  assemblages,
  method) {
  # Residual RF
  if (method %in% c("residual_rf", "sequential_residual_rf")) {
    return(names(assemblages)[grep("diff", names(assemblages))])
  }
  # Return for composition and abundance
  else if (predictor_variable %in% c("_composition","_abundance")) {
    return(names(assemblages)[grep("star", names(assemblages))])
  }
  # Return the feature in assemblage
  else if (predictor_variable %in% names(assemblages)) {
    return(predictor_variable)
  }
  # Return error otherwise
  else {
    print(paste("Error, predictor variable is not valid:", predictor_variable, method))
    return(NULL)
  }
}

convert_state_idx_to_vec <- function(
  num_species,
  state_index) {
  # Early exit if num_species > 32, since this method won't work
  if (num_species > 32) {
    stop("Number of species too large")
  }

  # Convert state_index binary mapping to rows with the state
  # state_index is 0-indexed
  return(as.integer(intToBits(state_index))[1:num_species])
}

convert_vec_to_state_idx <- function(
  num_species,
  state_vec) {
  # Convert binary state vector into integer (null state is 2^N)
  state_idx = 0
  state_existence = which(state_vec == 1)
  for (digits in state_existence) {
    state_idx = state_idx + 2^(digits-1)
  }
  if (state_idx == 0) {
    state_idx = 2^num_species
  }

  return(state_idx)
}

get_single_species_and_leave_one_out_state_idxs <- function(
  num_species) {
  # Numerically calculate index vectors
  single_species_state_idxs = 2^(0:(num_species-1))
  leave_one_out_state_idxs = 2^num_species - 1 - 2^(0:(num_species-1))
  
  return(union(single_species_state_idxs, leave_one_out_state_idxs))
}

get_full_state_grid <- function(
  num_species) {
  # Generate all states 
  full_states = do.call(
    rbind, 
    lapply(
      1:2^num_species, 
      function(x) convert_state_idx_to_vec(num_species, x)
    )
  )
  return(full_states)
}

get_state_assemblages_mapping <- function(
  num_species,
  assemblages) {
  # Return assemblage with the state index attached
  assemblages$state_idx = apply(
    assemblages[,1:num_species], 
    1, 
    function(x) convert_vec_to_state_idx(num_species, x)
  )

  return(assemblages)
}

generate_state_idxs_train <- function(
  experimental_design,
  num_train,
  assemblages,
  method,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Early exits for special cases
  if (method != "sequential_rf" && experimental_design == "prior") {
    return(NULL)
  }
  
  # Variable setup
  training_sample_size = num_train
  if (method == "sequential_rf") {
    training_sample_size = ceiling(num_train / (hyperparams$sequential_rf$iterations + 1))
  }
  # Extract the full states and the states that actually exist
  full_states = get_full_state_grid(num_species)
  existing_state_idxs = unique(assemblages[,'state_idx'])
  
  # Biased knowledge sampling from SP & LOO
  if (experimental_design == "prior") {
    prior_informed_training_idxs = get_single_species_and_leave_one_out_state_idxs(num_species)
    if (length(prior_informed_training_idxs) > num_train) {
      return(NULL)
    }
    # Early return is OK since this will always be less than the max size
    if (training_sample_size > length(prior_informed_training_idxs)) {
      sampling_idxs = setdiff(existing_state_idxs, prior_informed_training_idxs)
      rest_sampled_idxs = sample(
        x = sampling_idxs, 
        size = training_sample_size - length(prior_informed_training_idxs))
      return(union(prior_informed_training_idxs, rest_sampled_idxs))
    }
    else {
      return(prior_informed_training_idxs)
    }
  }
  # Uniform random sampling
  else if (experimental_design == "mixed") {
    sampling_idxs = existing_state_idxs 
  }
  # Sample from low richness states
  else if (grepl("low", experimental_design, fixed = TRUE)) {
    max_richness = as.numeric(gsub(
      "low-", "", as.character(experimental_design), fixed = TRUE))
    sampling_idxs = intersect(
      which(rowSums(full_states) <= max_richness),
      existing_state_idxs)
  }
  # Sample from high richness states
  else if (grepl("high", experimental_design, fixed = TRUE)) {
    min_richness = num_species - as.numeric(gsub(
      "high-", "" , as.character(experimental_design), fixed = TRUE))
    sampling_idxs = intersect(
      which(rowSums(full_states) >= min_richness),
      existing_state_idxs)
  }
  # Catch exception
  else {
    print(paste("Invalid training row sampling scheme:", experimental_design))
    return(NULL)
  }

  # Exit if requested training size is 0
  if (length(sampling_idxs) == 0) {
    print(paste(
      'No points for sampling, skipping schema:', 
      experimental_design
    ))
    return(NULL)
  }

  # Exit if requested training size is too large
  if (training_sample_size > length(sampling_idxs)) {
    print(paste(
      'Too many points requested for sampling, skipping schema:', 
      experimental_design
    ))
    return(NULL)
  }

  # Return sampled rows
  return(sample(x = sampling_idxs, size = training_sample_size))
}

generate_state_idxs_test <- function(
  experimental_design, 
  num_test,
  assemblages,
  num_species,
  skip_specific_state_idxs = NULL) {
  # The full list of testing candidates
  test_candidates_idxs = unique(assemblages[,'state_idx'])

  # If skipping specific states
  if (!is.null(skip_specific_state_idxs)) {
    # Take the set difference
    test_candidates_idxs = setdiff(
      test_candidates_idxs, skip_specific_state_idxs)
  }

  # Early break for edge cases
  num_filtered_test = min(num_test, length(test_candidates_idxs), 2^num_species)
  if (num_filtered_test == 0) {
    print(paste(
      'No test points generated, skipping schema:', 
      experimental_design
    ))
    return(NULL)
  }

  # Get the testing rows
  state_idxs_test = sample(
    x = test_candidates_idxs, 
    size = num_filtered_test)

  return(state_idxs_test)
}

get_rows_from_state_idxs <- function(
  state_idxs,
  assemblages) {
  # Filter the assemblages rows with state idx match
  return(which(assemblages[,'state_idx'] %in% state_idxs))
}

get_assemblages_subset_from_state_idxs <- function(
  state_idxs,
  assemblages) {
  # Filter the assemblages with state idx match
  return(assemblages %>% filter(state_idx %in% state_idxs))
}

process_data_features <- function(
  predictor_variable, 
  assemblages, 
  method,
  num_factor_bins = 9) {
  # Get columns for composition and abundance
  predictor_columns = get_predictor_columns(predictor_variable, assemblages, method)

  if (predictor_variable=="_composition") {
    # convert to presence/absence data - binary output
    assemblages[,predictor_columns] = (assemblages[,predictor_columns] > 0)
    for (col in predictor_columns) {
      assemblages[,col] = factor(assemblages[,col], levels = c(FALSE, TRUE))
    }
  }

  else if (predictor_variable=="_abundance") {
    # If we have abundance, convert the abundances to bins
    # assume ten total classes
    numeric_values = assemblages[,predictor_columns] %>% as.matrix %>% as.numeric
    quantile_bins = quantile(
      numeric_values[numeric_values > 0], 
      seq(0, 1, length.out = num_factor_bins), 
      na.rm = TRUE) 
    max_value = max(numeric_values, na.rm=TRUE)
    if (max_value == 0) {
      max_value = 1e-16 # a hack to get the breaks to work below
    }
    
    assemblages[,predictor_columns] = assemblages[,predictor_columns] %>% 
      mutate(across(everything(), function(x) {
        breaks = na.omit(as.numeric(unique(c(0, quantile_bins, max_value))))
        bin_means = (head(breaks, -1) + tail(breaks, -1))/2
        bin_means[1] = 0
        #print(bin_means) ### DEBUG
        return(cut(x, breaks=breaks,right=TRUE,include.lowest=TRUE,labels=bin_means))
      }
      ))
  }
  
  return(assemblages)
}

fit_rf_classifier_singlevar <- function(
  predictor_variable,
  assemblages,
  training_state_idxs,
  method,
  num_species) {
  # Get the training data
  data_training = get_assemblages_subset_from_state_idxs(
    training_state_idxs, assemblages)

  # Get the RF formula for fitting
  formula_rf_model = formula(sprintf(
    "%s~%s",
    predictor_variable,
    paste(names(data_training)[1:num_species], collapse = "+")
  ))

  # Get the case weights
  case_weights = freq_weight(as.numeric(data_training[,predictor_variable]))
  
  # Generate RF model
  rf_model = ranger(
    data = data_training,
    formula = formula_rf_model,
    importance = 'permutation',
    case.weights = case_weights,
    verbose = TRUE,
    num.trees = 500,
    mtry = ceiling(sqrt(num_species)),
    min.node.size = ceiling(sqrt(num_species))
  )

  rf_wrapper = list(
    "model" = rf_model,
    "training_state_idxs" = training_state_idxs
  )

  return(rf_wrapper)
}

fit_rf_classifier_multivar <- function(
  predictor_variable,
  assemblages,
  training_state_idxs,
  method,
  num_species) {
  # Get the training data
  data_training = get_assemblages_subset_from_state_idxs(
    training_state_idxs, assemblages)

  # Get columns for composition and abundance
  predictor_columns = get_predictor_columns(predictor_variable, assemblages, method)

  # Get the RF formula for fitting
  formula_rf_model = formula(sprintf(
    "Multivar(%s)~%s", 
    paste(predictor_columns, collapse=", "), 
    paste(names(data_training)[1:num_species], collapse="+")
  ))

  # Process data for multivar classification
  data_training_processed = process_data_features(
    predictor_variable, data_training, method)
  
  # Generate RF model
  rf_model = rfsrc(
    formula = formula_rf_model,
    data = data_training_processed,
    forest = TRUE,
    importance = "none",
    num.trees = 500,
    mtry = ceiling(sqrt(num_species)),
    min.node.size = ceiling(sqrt(num_species))
  )

  rf_wrapper = list(
    "model" = rf_model,
    "training_state_idxs" = training_state_idxs
  )

  return(rf_wrapper)
}

fit_rf_classifier <- function(
  predictor_variable,
  assemblages,
  training_state_idxs,
  method,
  num_species) {
  # Determine single vs. multi variable output
  is_single_var = !(predictor_variable %in% c("_abundance", "_composition"))

  # Pipeline for getting single vs. multi variable RF
  if (is_single_var) {
    return(fit_rf_classifier_singlevar(
      predictor_variable, assemblages, training_state_idxs, method, num_species))
  }
  else {
    return(fit_rf_classifier_multivar(
      predictor_variable, assemblages, training_state_idxs, method, num_species))
  }
}

normalize_scores_sequential_rf <- function(scores) {
  # Get the max and min
  max_score = max(scores)
  min_score = min(scores)

  # Normalization
  normalized_scores = (scores - min_score) / (max_score - min_score)

  return(normalized_scores)
}

estimate_uncertainty_sequential_rf <- function(
  predictor_variable,
  assemblages,
  training_state_idxs,
  candidate_state_idxs,
  method,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Define the bootstrap variables
  predictions_full = list()
  num_bootstrap = hyperparams$sequential_rf$bootstrap

  # Loop over each bootstrap fitting to obtain predictions
  for (iteration in 1:num_bootstrap) {
    # Bootstrap sample
    bootstrap_training_state_idxs = sample(
      training_state_idxs, replace = TRUE)
    
    # Fit rf model and predict
    bootstrap_rf_model = fit_rf_classifier(
      predictor_variable, assemblages, bootstrap_training_state_idxs, method, num_species)
    bootstrap_predictions = predict_rf_classifier(
      predictor_variable, bootstrap_rf_model, assemblages, candidate_state_idxs, num_species)

    # Append to boostrap predictions
    predictions_full[[iteration]] = bootstrap_predictions
    prediction_rows = nrow(bootstrap_predictions)
    prediction_dims = ncol(bootstrap_predictions)
  }

  # Get the uncertainty score - MSE of bootstrap prediction
  predictions_full_multidim = array(
    unlist(predictions_full), 
    dim = c(prediction_rows, prediction_dims, num_bootstrap))
  uncertainty_score = rowSums(apply(predictions_full_multidim, c(1, 2), var))
  
  # Normalize if needs be
  if (hyperparams$sequential_rf$normalize) {
    uncertainty_score = normalize_scores_sequential_rf(uncertainty_score) 
  }

  return(uncertainty_score)
}

estimate_diversity_sequential_rf <- function(
  candidate_state_idxs,
  assemblages, 
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Get the full state grid and extract states
  candidate_state = get_assemblages_subset_from_state_idxs(
    candidate_state_idxs, assemblages)[,1:num_species]
  
  # Get the nearest neighbors (k + 1 since self is included for cross-distance)
  nearest_k = hyperparams$sequential_rf$nearest_k
  nearest_neighbors_object = nn2(
    candidate_state, 
    query = candidate_state, 
    k = nearest_k + 1)
  
  # Get the diversity score (self-distance is 0 so we use k)
  diversity_score = rowSums(nearest_neighbors_object$nn.dists) / nearest_k

  # Normalize if needs be
  if (hyperparams$sequential_rf$normalize) {
    diversity_score = normalize_scores_sequential_rf(diversity_score) 
  }

  return(diversity_score)
}

estimate_density_sequential_rf <- function(
  training_state_idxs,
  candidate_state_idxs,
  assemblages, 
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Get the full state grid and extract states
  candidate_state = get_assemblages_subset_from_state_idxs(
    candidate_state_idxs, assemblages)[,1:num_species]
  training_state = get_assemblages_subset_from_state_idxs(
    training_state_idxs, assemblages)[,1:num_species]
  
  # Get the nearest neighbor
  nearest_neighbors_object = nn2(
    training_state, 
    query = candidate_state, 
    k = 1)
  
  # Get the density score
  density_score = nearest_neighbors_object$nn.dists[,1]

  # Normalize if needs be
  if (hyperparams$sequential_rf$normalize) {
    density_score = normalize_scores_sequential_rf(density_score) 
  }

  return(density_score)
}

estimate_best_candidates_sequential_rf <- function(
  predictor_variable,
  assemblages,
  training_state_idxs,
  batch_size, 
  method,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Set up comparison variables
  highest_score_min_heap = fibonacci_heap("numeric")
  candidate_state_idxs = setdiff(
      unique(assemblages[,'state_idx']), training_state_idxs)
  candidate_state = get_assemblages_subset_from_state_idxs(
    candidate_state_idxs, assemblages)[,1:num_species]
  score_weights = hyperparams$sequential_rf$score_weights
  
  # Get the scores
  uncertainty_score = estimate_uncertainty_sequential_rf(
    predictor_variable, assemblages, training_state_idxs, 
    candidate_state_idxs, method, num_species, hyperparams)
  diversity_score = estimate_diversity_sequential_rf(
    candidate_state_idxs, assemblages, num_species, hyperparams)
  density_score = estimate_density_sequential_rf(
    training_state_idxs, candidate_state_idxs, assemblages, 
    num_species, hyperparams)
  
  # Calculate the weighted score sum
  full_score = (
    uncertainty_score * score_weights$uncertainty +
    diversity_score * score_weights$diversity +
    density_score * score_weights$density
  )
  
  # Get the batch_size amount of best elements
  full_score_sorted = order(full_score, decreasing = TRUE)
  scored_state = candidate_state[full_score_sorted,]
  scored_state_with_idxs = get_state_assemblages_mapping(num_species, scored_state)
  scored_state_idxs = unique(scored_state_with_idxs[,'state_idx'])
  best_state_idxs = scored_state_idxs[1:min(batch_size, length(scored_state_idxs))]

  return(best_state_idxs)
}

get_batch_sizes <- function(
  num_samples,
  num_parts) {
  # Calculate the batch split
  lower_batch = floor(num_samples / num_parts)
  batch_remainder = (num_samples %% num_parts)
  batch_size_vec = rep(lower_batch, num_parts)
  if (batch_remainder > 0) {
    batch_size_vec[1:batch_remainder] = lower_batch + 1
  }
  return(batch_size_vec)
}

fit_sequential_rf_classifier <- function(
  predictor_variable,
  method,
  assemblages,
  num_train,
  training_state_idxs,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Initial setup of training rows and other variables
  sequential_training_state_idxs = training_state_idxs
  sequential_rf_iterations = hyperparams$sequential_rf$iterations
  max_states = length(unique(assemblages[,'state_idx']))
  batch_size_vec = get_batch_sizes(
    (num_train - length(training_state_idxs)),
    sequential_rf_iterations)

  for (iteration in 1:sequential_rf_iterations) {
    # Get the best improvement points
    batch_size = batch_size_vec[iteration]
    best_state_idxs = estimate_best_candidates_sequential_rf(
      predictor_variable, assemblages, 
      sequential_training_state_idxs, batch_size, method,
      num_species, hyperparams)
      
    # Update the sequential training rows
    sequential_training_state_idxs = union(
      sequential_training_state_idxs, best_state_idxs)
    if (length(sequential_training_state_idxs) == max_states) {break}
  }

  # Final fitting of random forest
  sequential_rf_wrapper = fit_rf_classifier(
    predictor_variable, assemblages, sequential_training_state_idxs, method, num_species)

  return(sequential_rf_wrapper)
}

vectorize_state_for_glv <- function (
  state_vec,
  initial_idx,
  num_species) {  
  # Early exit
  if (initial_idx > num_species || initial_idx < 1) {
    stop("Invalid species index for vectorization")
  }

  # Vectorize state for 1_n x N
  vectorized_state = rep(0, num_species^2)
  start_idx = (initial_idx - 1) * num_species + 1
  end_idx = initial_idx * num_species
  vectorized_state[start_idx:end_idx] = state_vec

  return(vectorized_state)
}

convert_training_row_to_glv_fitting <- function (
  assemblages_row,
  num_species) {
  # Set up the column labeling
  initial_cols = 1:num_species
  existence_cols = names(assemblages_row)[grep("star", names(assemblages_row))]
  initial_vec = as.numeric(assemblages_row[initial_cols])
  state_vec = as.numeric(assemblages_row[existence_cols])
  
  # Get the initial existence indices
  initial_existence = which(initial_vec == 1)

  # Vectorized state conversion
  glv_fitting_row = do.call(
    rbind, 
    lapply(
      initial_existence, 
      function(x) vectorize_state_for_glv(state_vec, x, num_species)
    )
  )

  return(glv_fitting_row)
}

create_glv_dataset <- function(
  assemblages,
  num_species) {
  # Get the vectorized versions
  glv_fitting = do.call(
    rbind,
    apply (
      assemblages,
      1,
      function(x) convert_training_row_to_glv_fitting(x, num_species)
    )
  )
  species_names = names(assemblages[1:num_species])
  colnames(glv_fitting) = crossing(species_names, species_names) %>% 
    unite("cols", sep=".") %>% 
    pull
  
  return(glv_fitting)
}

fit_glv_baseline <- function(
  assemblages,
  training_state_idxs,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Get the training data
  assemblages_training = get_assemblages_subset_from_state_idxs(
    training_state_idxs, assemblages)
  glv_training = create_glv_dataset(
    assemblages_training, num_species)

  # Fit a linear regression model
  # Basically, GLV steady state equation for each index boils down to
  # 1 + vec(A_matrix) * vec(1_n x N_h) = 0
  # so every time, we are trying to get the vector product to be close/equal to -1
  X = data.matrix(glv_training)
  y = data.matrix(rep(-1, nrow(glv_training)))
  A_effective_vec = as.numeric(ginv(X) %*% y)
  A_effective = matrix(A_effective_vec, nrow = num_species, byrow=TRUE)

  # Return the wrapper
  glv_wrapper = list(
    "model" = A_effective,
    "training_state_idxs" = training_state_idxs
  )

  return(glv_wrapper)
}

# fit_rf_regressor <- function() {

# }

fit_model <- function(
  predictor_variable,
  assemblages,
  num_train,
  training_state_idxs,
  num_species,
  method,
  hyperparams = MODEL_HYPERPARAMS) {
  # Switch for methods
  if (method == "naive") {
    # Naive fitting doesn't require model
    naive_wrapper = list(
      "model" = NULL,
      "training_state_idxs" = training_state_idxs
    )
    return(naive_wrapper)
  }
  else if (method == "rf") {
    # Random forest prediction
    return(fit_rf_classifier(
      predictor_variable, assemblages, training_state_idxs, method, num_species))
  }
  else if (method == "residual_rf") {
    # Random forest prediction
    return(fit_rf_regressor(
      predictor_variable, assemblages, training_state_idxs, method, num_species))
  }
  else if (method == "sequential_rf") {
    # Sequential fitting random forest prediction
    return(fit_sequential_rf_classifier(
      predictor_variable, method, assemblages, num_train, training_state_idxs, num_species))
  }
  else {
    print(paste("Invalid fitting method:", method))
    return(NULL)
  }
}

predict_rf_classifier_singlevar <- function(
  predictor_variable,
  rf_model,
  assemblages,
  predict_state_idxs,
  num_species) {
  # Get the prediction data
  data_predict = get_assemblages_subset_from_state_idxs(
    predict_state_idxs, assemblages)

  # Predict all values
  values_predicted = predict(
    rf_model$model, 
    data_predict
  )$predictions
  
  return(data.frame("pred"=values_predicted))
}

predict_rf_classifier_multivar <- function(
  predictor_variable,
  rf_model,
  assemblages,
  predict_state_idxs,
  num_species) {
  # Get the prediction data
  data_predict = get_assemblages_subset_from_state_idxs(
    predict_state_idxs, assemblages)[, 1:num_species]

  # Predict all values
  values_predicted_raw = predict(
    object = rf_model$model, 
    newdata = data_predict
  )
  # Process the predicted values
  if (predictor_variable == "_composition") {
    values_predicted = as.data.frame(sapply(
      values_predicted_raw$classOutput, 
      function(x) {as.logical(x$class)}, 
      simplify = FALSE
    ))
  }
  else if (predictor_variable == "_abundance") {
    # convert the class predictions back to numeric values
    values_predicted = as.data.frame(sapply(
      values_predicted_raw$classOutput, 
      function(x) {x$class}, 
      simplify = FALSE
    )) %>% 
      mutate(across(everything(), function(x) {as.numeric(as.character(x))}))
  }
  else {
    print("Error, predictor variable is not valid")
    return(NULL)
  }

  return(values_predicted)
}

predict_rf_classifier <- function(
  predictor_variable,
  rf_model,
  assemblages,
  predict_state_idxs,
  num_species) {
  # Determine single vs. multi variable output
  is_single_var = !(predictor_variable %in% c("_abundance", "_composition"))

  # Pipeline for getting single vs. multi variable RF
  if (is_single_var) {
    return(predict_rf_classifier_singlevar(
      predictor_variable, rf_model, assemblages, predict_state_idxs, num_species))
  }
  else {
    return(predict_rf_classifier_multivar(
      predictor_variable, rf_model, assemblages, predict_state_idxs, num_species))
  }
}

predict_naive_singlevar <- function(
  predictor_variable,
  assemblages,
  predict_state_idxs,
  num_species) {
  # Get the prediction data
  data_predict = get_assemblages_subset_from_state_idxs(
    predict_state_idxs, assemblages)

  # Predict all values
  if (is.factor(data_predict[, predictor_variable])) {
    # Randomly sample values of the factor from the training data
    values_predicted = sample(
      x = assemblages[,predictor_variable], 
      size = nrow(data_predict),
      replace = TRUE
    )
  }
  else {
    # Assume all species coexist
    if (predictor_variable == "richness") {
      # Use training data richness
      values_predicted = apply(
        data_predict[, 1:num_species], 1, sum)
    }
    # Pick mean value of the continuous variable
    else {
      values_predicted = mean(data_predict[, predictor_variable])
    }
  }
  
  return(values_predicted)
}

predict_naive_multivar <- function(
  predictor_variable,
  assemblages,
  predict_state_idxs,
  method,
  num_species) {
  # Get the prediction data
  data_predict = get_assemblages_subset_from_state_idxs(
    predict_state_idxs, assemblages)

  # Get columns for composition and abundance
  predictor_columns = get_predictor_columns(predictor_variable, assemblages, method)

  # If composition use the input presence/absences
  if (predictor_variable == "_composition") {
    values_predicted = (data_predict[, 1:num_species] > 0)
  }
  # If abundance use the mean training values masked by the input presence/absences
  else if (predictor_variable == "_abundance") {
    values_predicted = (data_predict[, 1:num_species] * 
      colMeans(data_predict[, predictor_columns]))
  }
  else {
    print(paste("Error, predictor variable is not valid:", predictor_variable))
    return(NULL)
  }

  return(values_predicted)
}

predict_naive <- function(
  predictor_variable,
  assemblages,
  predict_state_idxs,
  method,
  num_species) {
  # Determine single vs. multi variable output
  is_single_var = !(predictor_variable %in% c("_abundance", "_composition"))

  # Pipeline for getting single vs. multi variable RF
  if (is_single_var) {
    return(predict_naive_singlevar(
      predictor_variable, assemblages, predict_state_idxs, num_species))
  }
  else {
    return(predict_naive_multivar(
      predictor_variable, assemblages, predict_state_idxs, method, num_species))
  }
}

predict_glv_row <- function(
  glv_wrapper,
  assemblages,
  predict_idx,
  num_species) {
  # Get the prediction row and related quantities
  predict_row = assemblages[predict_idx,]
  predict_initial = as.numeric(predict_row[,1:num_species])
  initial_idx = which(predict_initial == 1)

  # Early exit if initial_idx has no initial state - GLV will default to 0
  if (length(initial_idx) == 0) {
    return(rep(0, num_species))
  }

  # Get the submatrix for the species that exist in the beginning
  A_effective_submatrix = glv_wrapper$model[initial_idx, initial_idx]

  # Make the prediction by -A^-1r, where r = 1 for us
  predict_existence_raw = as.numeric(
    -ginv(A_effective_submatrix) %*% rep(1, length(initial_idx)))
  
  # Populate vector by picking which species existed in the beginning
  predict_existence = rep(0, num_species)
  predict_existence[initial_idx] = clamp(predict_existence_raw, lower = 0)

  return(predict_existence)
}

predict_glv <- function(
  glv_wrapper,
  assemblages,
  predict_state_idxs,
  num_species) {
  # Get the prediction data
  data_predict = get_assemblages_subset_from_state_idxs(
    predict_state_idxs, assemblages)

  # Predict using fitted GLV model for each row
  values_predicted = do.call(
    rbind,
    lapply(
      1:nrow(data_predict),
      function(x) predict_glv_row(
        glv_wrapper, data_predict, x, num_species)
    )
  )

  return(values_predicted)
}

predict_model <- function(
  predictor_variable,
  rf_model,
  assemblages,
  predict_state_idxs,
  method,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Switch for methods
  if (method == "naive") {
    # Naive prediction
    return(predict_naive(
      predictor_variable, assemblages, predict_state_idxs, method, num_species))
  }
  else if (method == "rf" || method == "sequential_rf") {
    # Random forest prediction - same for sequential.
    # Since they both output same RF classifiers, just trained differently
    return(predict_rf_classifier(
      predictor_variable, rf_model, assemblages, predict_state_idxs, num_species))
  }
  else {
    print(paste("Invalid predicting method:", method))
    return(NULL)
  }
}

get_ground_truth_values <- function(
  predictor_variable,
  assemblages,
  predict_state_idxs,
  method,
  num_species) {
  # Get the prediction data
  data_predict = get_assemblages_subset_from_state_idxs(
    predict_state_idxs, assemblages)

  # Get columns for composition and abundance
  predictor_columns = get_predictor_columns(predictor_variable, assemblages, method)
  
  # Process the ground truth values
  if (predictor_variable == "_composition") {
    values_ground_truth = (data_predict[, 1:num_species] > 0)
  }
  else if (predictor_variable %in% c("_abundance", "richness")) {
    values_ground_truth = data_predict[, predictor_columns]
  }
  else {
    print(paste("Error, predictor variable is not valid:", predictor_variable))
    return(NULL)
  }

  return(values_ground_truth)
}

evaluate_balanced_accuracy_singlevar <- function(
  confusion_matrix) {
  # Directly use confusion matrix to evaluate balanced accuracy for single var
  return(confusion_matrix$byClass["Balanced Accuracy"])
}

evaluate_balanced_accuracy_multivar <- function(
  values_predicted,
  values_ground_truth) {
  # Get balanced accuracy with casewise mean
  balanced_accuracy = NA
  try(balanced_accuracy <- balanced_accuracy_casewise_mean(
    pred = values_predicted, 
    obs = values_ground_truth)
  )

  return(balanced_accuracy)
}

evaluate_mean_absolute_error_singlevar <- function(
  values_predicted,
  values_ground_truth) {
  return(mean_absolute_error(
    pred = as.numeric(unlist(values_predicted)), 
    obs = values_ground_truth
  ))
}

evaluate_mean_absolute_error_multivar <- function(
  values_predicted,
  values_ground_truth) {
  # Get MAE with casewise mean
  mean_absolute_error = NA
  try(mean_absolute_error <- mean_absolute_error_casewise_mean(
    pred = values_predicted, 
    obs = values_ground_truth)
  )

  return(mean_absolute_error)
}

evaluate_confusion_matrix <- function(
  values_predicted,
  values_ground_truth) {
  return(confusionMatrix(
    factor(values_predicted, levels = c(FALSE, TRUE)), 
    factor(values_ground_truth, levels = c(FALSE, TRUE))
  ))
}

evaluate_statistics <- function(
  values_predicted,
  values_ground_truth,
  predictor_variable,
  assemblages) {
  # Prepare all variables
  mean_absolute_error = NA
  balanced_accuracy = NA
  confusion_matrix = NA

  # Determine variable settings
  is_single_var = !(predictor_variable %in% c("_abundance", "_composition"))
  is_factor = FALSE
  if (is_single_var && 
    (is.factor(assemblages[,predictor_variable]) || 
      is.logical(assemblages[,predictor_variable]))
  ) {
    is_factor = TRUE
  }

  # Calculate all applicable statistics
  if (is_single_var) {
    if (is_factor) {
      confusion_matrix = evaluate_confusion_matrix(
        values_predicted, values_ground_truth)
      balanced_accuracy = evaluate_balanced_accuracy_singlevar(confusion_matrix)
    }
    else {
      mean_absolute_error = evaluate_mean_absolute_error_singlevar(
        values_predicted, values_ground_truth)
    }
  }
  else {
    if (predictor_variable == "_composition") {
      balanced_accuracy = evaluate_balanced_accuracy_multivar(
        values_predicted, values_ground_truth)
    }
    else {
      mean_absolute_error = evaluate_mean_absolute_error_multivar(
        values_predicted, values_ground_truth)
    }
  }

  return(list(
    mean_absolute_error = mean_absolute_error,
    balanced_accuracy = balanced_accuracy,
    confusion_matrix = confusion_matrix
  ))
}

clean_input_data <- function(input_file) {
  # Get the input data and create feasible and stable variable
  data = input_file %>%
    mutate(feasible.and.stable = factor(stable & feasible, levels=c(FALSE, TRUE))) %>%
    select(-feasible, -stable)
  
  # Flag quantile outliers
  data = quantile_max_trim(data)
  
  # Remove missing cases that arose from the above
  which_rows_na = data %>% 
    select(contains("star")) %>% 
    rowSums %>%
    is.na %>%
    which
  print(sprintf("Removed %d problematic case rows",length(which_rows_na)))
  data = data[-which_rows_na,]

  return(data)
}

generate_sample_size_sequences <- function(
  min_points,
  max_points,
  num_grid_points,
  num_data_row) {
  # Make the sample size sequence  
  sample_size_seq_all = unique(round(log_seq(
    min_points, max_points, length.out = num_grid_points)))

  # Trim to only the sizes that are compatible with the dataset
  sample_size_seq_all = sample_size_seq_all[sample_size_seq_all <= num_data_row]

  return(sample_size_seq_all)
}

perform_prediction_experiment_single <- function(
  predictor_variable,
  assemblages, 
  num_train,
  state_idxs_train,
  state_idxs_test,
  method, 
  num_species) {
  # If there is no training data, just exit early
  if (length(state_idxs_train)==0) {
    return(NULL)
  }

  # Fit a model with the wrapper
  fitted_model = fit_model(
    predictor_variable, assemblages, num_train, state_idxs_train, num_species, method)

  # Make predictions on both training and testing data
  model_predictions_train = predict_model(
    predictor_variable, fitted_model, assemblages, fitted_model$training_state_idxs, method, num_species)
  model_predictions_test = predict_model(
    predictor_variable, fitted_model, assemblages, state_idxs_test, method, num_species)

  # Early exit if predictions are not valid
  if (is.null(model_predictions_train) || is.null(model_predictions_test)) {
    return(NULL)
  }

  # Get ground truth labels and values
  ground_truth_train = get_ground_truth_values(
    predictor_variable, assemblages, fitted_model$training_state_idxs, method, num_species)
  ground_truth_test = get_ground_truth_values(
    predictor_variable, assemblages, state_idxs_test, method, num_species)
  
  # Gather statistics
  prediction_statistics_train = evaluate_statistics(
    model_predictions_train, ground_truth_train, predictor_variable, assemblages)
  prediction_statistics_test = evaluate_statistics(
    model_predictions_test, ground_truth_test, predictor_variable, assemblages)

  # Return full list of results
  return(list(
    model = fitted_model,
    pred_train = model_predictions_train,
    obs_train = ground_truth_train,
    cm_train = prediction_statistics_train$confusion_matrix,
    ba_train = prediction_statistics_train$balanced_accuracy,
    mae_train = prediction_statistics_train$mean_absolute_error,
    pred_test = model_predictions_test,
    obs_test = ground_truth_test,
    cm_test = prediction_statistics_test$confusion_matrix,
    ba_test = prediction_statistics_test$balanced_accuracy,
    mae_test = prediction_statistics_test$mean_absolute_error
  ))
}

write_to_csv_file <- function(
  file_to_save,
  directory_string,
  dataset_name,
  index,
  method,
  replicate_index,
  num_train,
  experimental_design,
  response,
  output) {
  # Prepare file name string
  csv_string = paste(
    directory_string, "/", dataset_name, "/",
    "index=", index, "_",
    "method=", method, "_",
    "rep=", replicate_index, "_",
    "num_train=", num_train, "_",
    "experimental_design=", experimental_design, "_",
    "response=", response, "_",
    "output=", output, ".csv",
    sep = ""
  )

  # Attempt to save to csv
  try(write.csv(file_to_save, file=csv_string, row.names = FALSE))
}

evaluate_initial_final_difference <- function(
  data_train,
  assemblages,
  num_species) {
  # Get outcomes
  initial_conditions = data_train[, 1:num_species] %>% as.matrix
  final_abundances = data_train %>% 
    select(contains("star")) %>% as.matrix
  
  # Count # of species that were present but went absent
  num_losses_mean = mean(apply(
    (final_abundances==0) & (initial_conditions==1), 1, sum, na.rm = TRUE))

  # Figure out abundance distribution in training
  abundance_final_skewness_mean = skewness(
    as.numeric(final_abundances), na.rm = TRUE)
  abundance_final_skewness_nonzero_mean = skewness(
    as.numeric(final_abundances)[as.numeric(final_abundances) > 0], na.rm = TRUE)
  
  # Determine the overall dataset 95% abundance quantile
  abundances_dataset_all = assemblages %>% 
    select(contains("star")) %>% as.matrix %>% as.numeric
  abundance_q95_dataset = quantile(abundances_dataset_all, 0.95, na.rm = TRUE)

  # Determine overall dataset skewness
  abundance_skewness_dataset = skewness(abundances_dataset_all, na.rm = TRUE)
  abundance_skewness_nonzero_dataset = skewness(
    abundances_dataset_all[abundances_dataset_all > 0], na.rm = TRUE)

  return(list(
    num_losses_mean = num_losses_mean,
    abundance_final_skewness_mean = abundance_final_skewness_mean,
    abundance_final_skewness_nonzero_mean = abundance_final_skewness_nonzero_mean,
    abundance_q95_dataset = abundance_q95_dataset,
    abundance_skewness_dataset = abundance_skewness_dataset,
    abundance_skewness_nonzero_dataset = abundance_skewness_nonzero_dataset
  ))
}

perform_prediction_experiment_parallel_wrapper <- function(
  directory_string,
  dataset_name,
  num_species,
  num_replicates_in_data,
  index,
  assemblages,
  results_table) {
  # Extract variables
  num_train = results_table$num_train[index]
  num_test = results_table$num_test[index]
  replicate_index = results_table$replicate_index[index]
  method = results_table$method[index]
  experimental_design = results_table$experimental_design[index]
  existing_state_idxs = unique(assemblages[,'state_idx'])
  
  # Subsample assemblage with desired replicates
  assemblages = assemblages %>% 
    slice_sample(n = num_replicates_in_data, by = state_idx)

  # Get training & testing set
  state_idxs_train = generate_state_idxs_train( 
    experimental_design, num_train, assemblages, method, num_species)
  state_idxs_test = generate_state_idxs_test( 
    experimental_design, num_test, assemblages, num_species)
  results_table$num_test[index] = length(state_idxs_test)

  # Early exit
  if (is.null(state_idxs_train) || is.null(state_idxs_test)) {
    return(NULL)
  }
  else if (length(existing_state_idxs) < num_train) {
    return(NULL)
  }

  # Get the rows for printing out
  data_train = get_assemblages_subset_from_state_idxs(
    state_idxs_train, assemblages)
  data_test = get_assemblages_subset_from_state_idxs(
    state_idxs_test, assemblages)

  # Save train and test data to file
  write_to_csv_file(
    data_train, directory_string, dataset_name, 
    index, method, replicate_index, num_train, 
    experimental_design, "abundance", "experiment_train")
  write_to_csv_file(
    data_test, directory_string, dataset_name, 
    index, method, replicate_index, num_train, 
    experimental_design, "abundance", "experiment_test")
  
  # Perform experiments for abundance, composition, richness
  for (response in c("_abundance", "_composition", "richness")) {
    response_save = gsub("_", "", response, fixed = TRUE)
    experiment_result = perform_prediction_experiment_single(
      response, assemblages, num_train, 
      state_idxs_train, state_idxs_test,
      method, num_species)
    
    # Early exit if results are not valid
    if (is.null(experiment_result)) {
      print("Returning NULL result")
      return(NULL)
    }

    # Get the proper diagnostics
    if (response == "_abundance") {
      results_table$abundance_mae_mean_train[index] = experiment_result$mae_train
      results_table$abundance_mae_mean_test[index] = experiment_result$mae_test
    }
    else if (response == "_composition") {
      results_table$composition_balanced_accuracy_mean_train[index] = experiment_result$ba_train
      results_table$composition_balanced_accuracy_mean_test[index] = experiment_result$ba_test
    }
    else if (response == "richness") {
      results_table$richness_mae_train[index] = experiment_result$mae_train
      results_table$richness_mae_test[index] = experiment_result$mae_test
    }

    # Save relevant variables
    write_to_csv_file(
      experiment_result$pred_train, directory_string, dataset_name, 
      index, method, replicate_index, num_train, 
      experimental_design, response_save, "pred_train")
    write_to_csv_file(
      experiment_result$pred_test, directory_string, dataset_name, 
      index, method, replicate_index, num_train, 
      experimental_design, response_save, "pred_test")
    write_to_csv_file(
      experiment_result$obs_train, directory_string, dataset_name, 
      index, method, replicate_index, num_train, 
      experimental_design, response_save, "obs_train")
    write_to_csv_file(
      experiment_result$obs_test, directory_string, dataset_name, 
      index, method, replicate_index, num_train, 
      experimental_design, response_save, "obs_test")
    
    # Save the final training set if sequential
    if (method == "sequential_rf") {
      data_train_sequential = get_assemblages_subset_from_state_idxs(
        experiment_result$model$training_state_idxs, assemblages)
      write_to_csv_file(
        data_train_sequential, 
        directory_string, dataset_name, 
        index, method, replicate_index, num_train, 
        experimental_design, "abundance", 
        paste("experiment_sequential_train", response, sep = ""))
    }
  }

  # Get the difference statistics
  difference_stats_train = evaluate_initial_final_difference(
    data_train, assemblages, num_species)
  
  # Record to table
  results_table$num_losses_mean[index] = difference_stats_train$num_losses_mean
  results_table$abundance_final_skewness_mean[index] = difference_stats_train$abundance_final_skewness_mean
  results_table$abundance_final_skewness_nonzero_mean[index] = difference_stats_train$abundance_final_skewness_nonzero_mean
  results_table$abundance_q95_dataset[index] = difference_stats_train$abundance_q95_dataset
  results_table$abundance_skewness_dataset[index] = difference_stats_train$abundance_skewness_dataset
  results_table$abundance_skewness_nonzero_dataset[index] = difference_stats_train$abundance_skewness_nonzero_dataset

  return(results_table[index, , drop = FALSE])
}

perform_prediction_experiment_full <- function(
  directory_string,
  input_file,
  dataset_name,
  num_species, 
  method_list,
  experimental_design_list,
  num_replicates_in_data = 1, 
  num_test = NUM_TEST,
  num_replicates_in_fitting = REPLICATES,
  num_grid_points = GRID_POINTS,
  min_points = MIN_POINTS,
  max_points = MAX_POINTS,
  parallelized = (CORES > 1)) {
  # Clean data and prep
  assemblages = clean_input_data(input_file)
  assemblages = get_state_assemblages_mapping(num_species, assemblages)
  sample_size_seq_all = generate_sample_size_sequences(
    min_points, max_points, num_grid_points, nrow(assemblages))
  dir.create(file.path(directory_string, dataset_name), recursive = TRUE)

  # Make the giant results table
  results_table = expand.grid(replicate_index=1:num_replicates_in_fitting, 
                              method=method_list,
                              experimental_design=experimental_design_list,
                              num_train=sample_size_seq_all,
                              num_test=num_test, 
                              richness_mae_test=NA,
                              feasible_and_stable_balanced_accuracy_test=NA,
                              composition_balanced_accuracy_mean_test=NA,
                              abundance_mae_mean_test=NA,
                              richness_mae_train=NA,
                              feasible_and_stable_balanced_accuracy_train=NA,
                              composition_balanced_accuracy_mean_train=NA,
                              abundance_mae_mean_train=NA,
                              num_species_dataset=NA,
                              num_replicates_dataset=NA,
                              num_cases_dataset=NA,
                              num_losses_mean_train=NA,
                              abundance_q95_dataset=NA,
                              abundance_skewness_dataset=NA,
                              abundance_skewness_nonzero_dataset=NA,
                              abundance_final_skewness_mean_train=NA,
                              abundance_final_skewness_nonzero_mean_train=NA
  )

  # Apply multi core parallelization
  indices = 1:nrow(results_table)
  if (parallelized) {
    results_list = mclapply(indices, function(index) {
      perform_prediction_experiment_parallel_wrapper(
        directory_string, dataset_name, num_species, 
        num_replicates_in_data, index, 
        assemblages, results_table)
    }, mc.cores = CORES, mc.preschedule = FALSE)
  }
  else {
    print(length(indices))
    for (index in indices) {
      print(index)
      print(results_table[index,])
      perform_prediction_experiment_parallel_wrapper(
        directory_string, dataset_name, num_species, 
        num_replicates_in_data, index, 
        assemblages, results_table)
    }
    # results_list = lapply(indices, function(index) {
    #   perform_prediction_experiment_parallel_wrapper(
    #     directory_string, dataset_name, num_species, 
    #     num_replicates_in_data, index, 
    #     assemblages, results_table)
    # })
  }
  # Get indices with errors
  indices_errors = which(sapply(results_list, class) == "try-error")
  indices_skipped = which(sapply(results_list, class) == "null")
  indices_good = which(sapply(results_list, class) == "data.frame")
  
  # Print indices with errors
  if (length(indices_errors) > 0) {
    print(paste("Indices with errors: ", paste(indices_errors, collapse = ', ')))
    print(results_list[indices_errors])
  }

  # Write results table
  results_df = NULL
  try(results_df <- rbindlist(results_list[indices_good]))
  if (!is.null(results_df)) {
    write.csv(
      results_df, file=sprintf('outputs/statistical/results_%s.csv',dataset_name), 
      row.names=FALSE)
  }

  # Save the raw output too in case of a rbind issue for error checking
  saveRDS(results_list, file = sprintf(
    'outputs/statistical/results_%s.Rdata', dataset_name))
  
  return(results_list)  
}

