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

# Source utility functions
source('utils/freq_weight.R')
source('utils/skill_statistics.R')
source('utils/log_seq.R')
source('utils/quantile_trim.R')


get_predictor_columns <- function(
  predictor_variable,
  assemblages) {
  # Return for composition and abundance
  if (predictor_variable %in% c("_composition","_abundance")) {
    return(names(assemblages)[grep("star", names(assemblages))])
  }
  # Return the feature in assemblage
  else if (predictor_variable %in% names(assemblages)) {
    return(predictor_variable)
  }
  # Return error otherwise
  else {
    print(paste("Error, predictor variable is not valid:", predictor_variable))
    return(NULL)
  }
}

get_single_species_and_leave_one_out_rows <- function(
  assemblages, 
  num_species) {
  # Get the single species only rows and leave one out rows
  single_species_rows = which(rowSums(assemblages[,1:num_species]) == 1)
  leave_one_out_rows = which(rowSums(assemblages[,1:num_species]) == num_species - 1)

  return(union(single_species_rows, leave_one_out_rows))
}

generate_rows_train <- function(
  experimental_design,
  method,
  num_train,
  assemblages,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Early exits for special cases
  if (method != "sequential_rf" && experimental_design == "prior") {
    return(NULL)
  }
  
  # Variable setup
  initial_richness = apply(assemblages[,1:num_species], 1, sum)
  max_rows = nrow(assemblages)
  training_sample_size = num_train
  if (method == "sequential_rf") {
    training_sample_size = floor(num_train / (hyperparams$sequential_rf$iterations + 1))
  }
  
  # Biased knowledge sampling from SP & LOO
  if (experimental_design == "prior") {
    prior_informed_training_rows = get_single_species_and_leave_one_out_rows(
      assemblages, num_species)
    if (length(prior_informed_training_rows) > num_train) {
      return(NULL)
    }
    # Early return is OK since this will always be less than the max size
    if (training_sample_size > length(prior_informed_training_rows)) {
      sampling_rows = setdiff(1:max_rows, prior_informed_training_rows)
      rest_sampled_rows = sample(
        x = sampling_rows, 
        size = training_sample_size - length(prior_informed_training_rows))
      return(union(prior_informed_training_rows, rest_sampled_rows))
    }
    else {
      return(prior_informed_training_rows)
    }
  }
  # Uniform random sampling
  else if (experimental_design == "mixed") {
    sampling_rows = 1:max_rows 
  }
  # Sample from low richness rows
  else if (grepl("low", experimental_design, fixed = TRUE)) {
    max_richness = as.numeric(gsub(
      "low-", "", as.character(experimental_design), fixed = TRUE))
    sampling_rows = which(initial_richness <= max_richness)
  }
  # Sample from high richness rows
  else if (grepl("high", experimental_design, fixed = TRUE)) {
    min_richness = num_species - as.numeric(gsub(
      "high-", "" , as.character(experimental_design), fixed = TRUE))
    sampling_rows = which(initial_richness >= min_richness)
  }
  # Catch exception
  else {
    print(paste("Invalid training row sampling scheme:", experimental_design))
    return(NULL)
  }

  # Exit if requested training size is too large
  if (training_sample_size > length(sampling_rows)) {
    print(paste(
      'Too many points requested for sampling, skipping schema:', 
      experimental_design
    ))
    return(NULL)
  }

  # Return sampled rows
  return(sample(x = sampling_rows, size = training_sample_size))
}

generate_rows_test <- function(
  experimental_design, 
  assemblages, 
  num_test,
  num_species,
  skip_specific_states_rows = NULL) {
  # The full list of testing candidates
  test_candidates_indices = 1:nrow(assemblages)

  # If skipping specific states (for instance, single species & leave one out rows)
  if (!is.null(skip_specific_states_rows)) {
    # Take the set difference
    test_candidates_indices = setdiff(
      test_candidates_indices, skip_specific_states_rows)
  }

  # Early break for edge cases
  num_filtered_test = min(num_test, length(test_candidates_indices), nrow(assemblages))
  if (num_filtered_test == 0) {
    print(paste(
      'No test points generated, skipping schema:', 
      experimental_design
    ))
    return(NULL)
  }

  # Get the testing rows
  rows_test = sample(
    x = test_candidates_indices, 
    size = num_filtered_test)

  return(rows_test)
}

process_data_features <- function(
  predictor_variable, 
  assemblages, 
  num_factor_bins = 9) {
  # Get columns for composition and abundance
  predictor_columns = get_predictor_columns(predictor_variable, assemblages)

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
  training_rows,
  num_species) {
  # Get training data
  data_training = assemblages[training_rows,]

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
    "training_rows" = training_rows
  )

  return(rf_wrapper)
}

fit_rf_classifier_multivar <- function(
  predictor_variable,
  assemblages,
  training_rows,
  num_species) {
  # Get training data
  data_training = assemblages[training_rows,]

  # Get columns for composition and abundance
  predictor_columns = get_predictor_columns(predictor_variable, data_training)

  # Get the RF formula for fitting
  formula_rf_model = formula(sprintf(
    "Multivar(%s)~%s", 
    paste(predictor_columns, collapse=", "), 
    paste(names(data_training)[1:num_species], collapse="+")
  ))

  # Process data for multivar classification
  data_training_processed = process_data_features(
    predictor_variable, data_training)
  
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
    "training_rows" = training_rows
  )

  return(rf_wrapper)
}

fit_rf_classifier <- function(
  predictor_variable,
  assemblages,
  training_rows,
  num_species) {
  # Determine single vs. multi variable output
  is_single_var = !(predictor_variable %in% c("_abundance", "_composition"))

  # Pipeline for getting single vs. multi variable RF
  if (is_single_var) {
    return(fit_rf_classifier_singlevar(
      predictor_variable, assemblages, training_rows, num_species))
  }
  else {
    return(fit_rf_classifier_multivar(
      predictor_variable, assemblages, training_rows, num_species))
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
  training_rows,
  candidate_rows,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Define the bootstrap variables
  predictions_full = list()
  num_bootstrap = hyperparams$sequential_rf$bootstrap

  # Loop over each bootstrap fitting to obtain predictions
  for (iteration in 1:num_bootstrap) {
    # Bootstrap sample
    bootstrap_training_rows = sample(
      training_rows, replace = TRUE)

    # Fit rf model and predict
    bootstrap_rf_model = fit_rf_classifier(
      predictor_variable, assemblages, bootstrap_training_rows, num_species)
    bootstrap_predictions = predict_rf_classifier(
      predictor_variable, bootstrap_rf_model, assemblages, candidate_rows, num_species)
  
    # Append to boostrap predictions
    predictions_full[[iteration]] = bootstrap_predictions
    prediction_dims = ncol(bootstrap_predictions)
  }

  # Get the uncertainty score - MSE of bootstrap prediction
  predictions_full_multidim = array(
    unlist(predictions_full), 
    dim = c(length(candidate_rows), prediction_dims, num_bootstrap))
  uncertainty_score = rowSums(apply(predictions_full_multidim, c(1, 2), var))
  
  # Normalize if needs be
  if (hyperparams$sequential_rf$normalize) {
    uncertainty_score = normalize_scores_sequential_rf(uncertainty_score) 
  }

  return(uncertainty_score)
}

estimate_diversity_sequential_rf <- function(
  select_candidate_data,
  hyperparams = MODEL_HYPERPARAMS) {
  # Get the nearest neighbors (k + 1 since self is included for cross-distance)
  nearest_k = hyperparams$sequential_rf$nearest_k
  nearest_neighbors_object = nn2(
    select_candidate_data, 
    query = select_candidate_data, 
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
  select_training_data,
  select_candidate_data,
  hyperparams = MODEL_HYPERPARAMS) {
  # Get the nearest neighbor
  nearest_neighbors_object = nn2(
    select_training_data, 
    query = select_candidate_data, 
    k = 1)
  
  # Get the density score (self-distance is 0 so we use k)
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
  training_rows,
  batch_size,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Set up comparison variables
  highest_score_min_heap = fibonacci_heap("numeric")
  candidate_rows = setdiff(
      1:nrow(assemblages), training_rows)
  # TODO: Need some way of getting candidates so that the assemblage doesn't overlap...
  score_weights = hyperparams$sequential_rf$score_weights
  select_training_data = assemblages[training_rows, 1:num_species]
  select_candidate_data = assemblages[candidate_rows, 1:num_species]

  # Get the scores
  uncertainty_score = estimate_uncertainty_sequential_rf(
    predictor_variable, assemblages, training_rows, 
    candidate_rows, num_species, hyperparams)
  diversity_score = estimate_diversity_sequential_rf(
    select_candidate_data, hyperparams)
  density_score = estimate_density_sequential_rf(
    select_training_data, select_candidate_data, hyperparams)
  
  # Calculate the weighted score sum
  full_score = (
    uncertainty_score * score_weights$uncertainty +
    diversity_score * score_weights$diversity +
    density_score * score_weights$density
  )

  # Get the batch_size amount of best elements (nice sorting trick)
  full_score_sorted = order(full_score, decreasing = FALSE)
  best_row_indices = candidate_rows[full_score_sorted[1:batch_size]]
  
  return(best_row_indices)
}

fit_sequential_rf_classifier <- function(
  predictor_variable,
  assemblages,
  num_train,
  training_rows,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Initial setup of training rows and other variables
  sequential_training_rows = training_rows
  sequential_rf_iterations = hyperparams$sequential_rf$iterations
  batch_size = floor(
    (num_train - length(training_rows)) / 
    sequential_rf_iterations
  )

  for (iteration in 1:sequential_rf_iterations) {
    # Get the best improvement points
    best_row_indices = estimate_best_candidates_sequential_rf(
      predictor_variable, assemblages, 
      sequential_training_rows, batch_size, num_species, hyperparams)
    # Update the sequential training rows
    sequential_training_rows = union(
      sequential_training_rows, best_row_indices)
  }

  # Final fitting of random forest
  sequential_rf_wrapper = fit_rf_classifier(
    predictor_variable, assemblages, sequential_training_rows, num_species)

  return(sequential_rf_wrapper)
}

fit_model <- function(
  predictor_variable,
  assemblages,
  num_train,
  training_rows,
  num_species,
  method,
  hyperparams = MODEL_HYPERPARAMS) {
  # Switch for methods
  if (method == "naive") {
    # Naive fitting doesn't require model
    naive_wrapper = list(
      "model" = NULL,
      "training_rows" = training_rows
    )
    return(naive_wrapper)
  }
  else if (method == "rf") {
    # Random forest prediction
    return(fit_rf_classifier(
      predictor_variable, assemblages, training_rows, num_species))
  }
  else if (method == "sequential_rf") {
    # Sequential fitting random forest prediction
    return(fit_sequential_rf_classifier(
      predictor_variable, assemblages, num_train, training_rows, num_species))
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
  predict_rows,
  num_species) {
  # Predict all values
  values_predicted = predict(
    rf_model$model, 
    assemblages[predict_rows,]
  )$predictions
  
  return(values_predicted)
}

predict_rf_classifier_multivar <- function(
  predictor_variable,
  rf_model,
  assemblages,
  predict_rows,
  num_species) {
  # Predict all values
  values_predicted_raw = predict(
    object = rf_model$model, 
    newdata = assemblages[predict_rows, 1:num_species]
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
  predict_rows,
  num_species) {
  # Determine single vs. multi variable output
  is_single_var = !(predictor_variable %in% c("_abundance", "_composition"))

  # Pipeline for getting single vs. multi variable RF
  if (is_single_var) {
    return(predict_rf_classifier_singlevar(
      predictor_variable, rf_model, assemblages, predict_rows, num_species))
  }
  else {
    return(predict_rf_classifier_multivar(
      predictor_variable, rf_model, assemblages, predict_rows, num_species))
  }
}

predict_naive_singlevar <- function(
  predictor_variable,
  assemblages,
  predict_rows,
  num_species) {
  # Predict all values
  if (is.factor(assemblages[predict_rows, predictor_variable])) {
    # Randomly sample values of the factor from the training data
    values_predicted = sample(
      x = assemblages[,predictor_variable], 
      size = length(predict_rows),
      replace = TRUE
    )
  }
  else {
    # Assume all species coexist
    if (predictor_variable == "richness") {
      # Use training data richness
      values_predicted = apply(
        assemblages[predict_rows, 1:num_species], 1, sum)
    }
    # Pick mean value of the continuous variable
    else {
      values_predicted = mean(assemblages[predict_rows, predictor_variable])
    }
  }
  
  return(values_predicted)
}

predict_naive_multivar <- function(
  predictor_variable,
  assemblages,
  predict_rows,
  num_species) {
  # Get columns for composition and abundance
  predictor_columns = get_predictor_columns(predictor_variable, assemblages)

  # If composition use the input presence/absences
  if (predictor_variable == "_composition") {
    values_predicted = (assemblages[predict_rows, 1:num_species] > 0)
  }
  # If abundance use the mean training values masked by the input presence/absences
  else if (predictor_variable == "_abundance") {
    values_predicted = (assemblages[predict_rows, 1:num_species] * 
      colMeans(assemblages[predict_rows, predictor_columns]))
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
  predict_rows,
  num_species) {
  # Determine single vs. multi variable output
  is_single_var = !(predictor_variable %in% c("_abundance", "_composition"))

  # Pipeline for getting single vs. multi variable RF
  if (is_single_var) {
    return(predict_naive_singlevar(
      predictor_variable, assemblages, predict_rows, num_species))
  }
  else {
    return(predict_naive_multivar(
      predictor_variable, assemblages, predict_rows, num_species))
  }
}

predict_model <- function(
  predictor_variable,
  rf_model,
  assemblages,
  predict_rows,
  method,
  num_species,
  hyperparams = MODEL_HYPERPARAMS) {
  # Switch for methods
  if (method == "naive") {
    # Naive prediction
    return(predict_naive(
      predictor_variable, assemblages, predict_rows, num_species))
  }
  else if (method == "rf" || method == "sequential_rf") {
    # Random forest prediction - same for sequential.
    # Since they both output same RF classifiers, just trained differently
    return(predict_rf_classifier(
      predictor_variable, rf_model, assemblages, predict_rows, num_species))
  }
  else {
    print(paste("Invalid predicting method:", method))
    return(NULL)
  }
}

get_ground_truth_values <- function(
  predictor_variable,
  assemblages,
  predict_rows,
  num_species) {
  # Get columns for composition and abundance
  predictor_columns = get_predictor_columns(predictor_variable, assemblages)
  
  # Process the ground truth values
  if (predictor_variable == "_composition") {
    values_ground_truth = (assemblages[predict_rows, 1:num_species] > 0)
  }
  else if (predictor_variable %in% c("_abundance", "richness")) {
    values_ground_truth = assemblages[predict_rows, predictor_columns]
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
    pred = values_predicted, 
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
  rows_train,
  rows_test, 
  method, 
  num_species) {
  # If there is no training data, just exit early
  if (length(rows_train)==0) {
    return(NULL)
  }

  # Get training & testing data
  data_train = assemblages[rows_train,]
  data_test = assemblages[rows_test,]

  # Fit a model with the wrapper
  fitted_model = fit_model(
    predictor_variable, assemblages, num_train, rows_train, num_species, method)

  # Make predictions on both training and testing data
  model_predictions_train = predict_model(
    predictor_variable, fitted_model, assemblages, fitted_model$training_rows, method, num_species)
  model_predictions_test = predict_model(
    predictor_variable, fitted_model, assemblages, rows_test, method, num_species)

  # Early exit if predictions are not valid
  if (is.null(model_predictions_train) || is.null(model_predictions_test)) {
    return(NULL)
  }

  # Get ground truth labels and values
  ground_truth_train = get_ground_truth_values(
    predictor_variable, assemblages, fitted_model$training_rows, num_species)
  ground_truth_test = get_ground_truth_values(
    predictor_variable, assemblages, rows_test, num_species)
  
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
  evaluate_rows,
  assemblages,
  num_species) {
  # Get outcomes
  initial_conditions = assemblages[evaluate_rows, 1:num_species] %>% as.matrix
  final_abundances = assemblages[evaluate_rows,] %>% 
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
  index,
  assemblages,
  results_table) {
  # Extract variables
  num_train = results_table$num_train[index]
  num_test = results_table$num_test[index]
  replicate_index = results_table$replicate_index[index]
  method = results_table$method[index]
  experimental_design = results_table$experimental_design[index]
  
  # Get training & testing set
  rows_train = generate_rows_train(
    experimental_design, method, num_train, assemblages, num_species)
  rows_test = generate_rows_test(
    experimental_design, assemblages, num_test, num_species,
    skip_specific_states_rows = rows_train)
  results_table$num_test[index] = length(rows_test)

  # Early exit
  if (is.null(rows_train) || is.null(rows_test)) {
    return(NULL)
  }

  # Save train and test data to file
  write_to_csv_file(
    assemblages[rows_train,], directory_string, dataset_name, 
    index, method, replicate_index, num_train, 
    experimental_design, "abundance", "experiment_train")
  write_to_csv_file(
    assemblages[rows_test,], directory_string, dataset_name, 
    index, method, replicate_index, num_train, 
    experimental_design, "abundance", "experiment_test")
  
  # Perform experiments for abundance, composition, richness
  for (response in c("_abundance", "_composition", "richness")) {
    response_save = gsub("_", "", response, fixed = TRUE)
    experiment_result = perform_prediction_experiment_single(
      response, assemblages, num_train, rows_train, rows_test, 
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
      write_to_csv_file(
        assemblages[experiment_result$model$training_rows,], 
        directory_string, dataset_name, 
        index, method, replicate_index, num_train, 
        experimental_design, "abundance", 
        paste("experiment_sequential_train", response, sep = ""))
    }
  }

  # Get the difference statistics
  difference_stats_train = evaluate_initial_final_difference(
    rows_train, assemblages, num_species)
  
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
  sample_size_seq_all = generate_sample_size_sequences(
    min_points, max_points, num_grid_points, nrow(assemblages))
  dir.create(file.path(directory_string, dataset_name), recursive = TRUE)

  # Make the giant results table
  results_table = expand.grid(replicate_index=1:num_replicates_in_fitting, 
                              method=method_list,
                              richness_mae_test=NA,
                              feasible_and_stable_balanced_accuracy_test=NA,
                              composition_balanced_accuracy_mean_test=NA,
                              abundance_mae_mean_test=NA,
                              richness_mae_train=NA,
                              feasible_and_stable_balanced_accuracy_train=NA,
                              composition_balanced_accuracy_mean_train=NA,
                              abundance_mae_mean_train=NA,
                              experimental_design=experimental_design_list,
                              num_train=sample_size_seq_all,
                              num_test=num_test, 
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
        index, assemblages, results_table)
    }, mc.cores = CORES)
  }
  else {
    results_list = lapply(indices, function(index) {
      perform_prediction_experiment_parallel_wrapper(
        directory_string, dataset_name, num_species, 
        index, assemblages, results_table)
    })
  }
  
  # Write results table
  results_df = NULL
  try(results_df <- rbindlist(results_list))
  if (!is.null(results_df)) {
    write.csv(
      results_df, file=sprintf('outputs/statistical/results_%s.csv',dataset_name), 
      row.names=FALSE)
  }

  # Save the raw output too in case of a rbind issue for error checking
  saveRDS(results_list, file = sprintf(
    'outputs/statistical/results_%s.Rdata', dataset_name))
  
  # Get indices with errors
  indices_errors = which(sapply(results_list, class) == "try-error")
  if (length(indices_errors) > 0) {
    print(data.frame(indices_errors, sapply(indices_errors, function(j) { 
      attr(results_list[[j]], "condition")
    })))
  }
  
  return(results_list)  
}

