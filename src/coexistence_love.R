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
    return(names(data_for_rf_training)[grep("star", names(data_for_rf_training))])
  }
  # Return the feature in assemblage
  else if (predictor_variable %in% names(assemblages)) {
    return(predictor_variable)
  }
  # Return error otherwise
  else {
    stop("Error, predictor variable is not valid")
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

generate_rows_test <- function(
  assemblages, 
  num_test,
  num_species,
  skip_specific_states = NULL) {
  # The full list of testing candidates
  test_candidates_indices = 1:nrow(assemblages)

  # If skipping specific states (for instance, single species & leave one out rows)
  if (!is.null(skip_specific_states)) {
    # Take the set difference
    test_candidates_indices = setdiff(test_candidates_indices, skip_specific_states)
  }

  # Get the testing rows
  rows_test = sample(x = test_candidates_indices, size = num_test)

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
      assemblages[,col] = factor(assemblages[,col],levels = c(FALSE, TRUE))
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

fit_rf_singlevar <- function(
  predictor_variable,
  data_training,
  num_species) {
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

  return(rf_model)
}

fit_rf_multivar <- function(
  predictor_variable,
  data_training,
  num_species) {
  # Get columns for composition and abundance
  predictor_columns = get_predictor_columns(predictor_variable, assemblages)

  # Get the RF formula for fitting
  formula_rf_model = formula(sprintf(
    "Multivar(%s)~%s", 
    paste(predictor_columns, collapse=", "), 
    paste(names(data_for_rf_training)[1:num_species], collapse="+")
  ))
  
  # Generate RF model
  rf_model = rfsrc(
    formula = formula_rf_model,
    data = data_for_rf_training,
    forest = TRUE,
    importance = "none",
    num.trees = 500,
    mtry = ceiling(sqrt(num_species)),
    min.node.size = ceiling(sqrt(num_species))
  )

  return(rf_model)
}

fit_model <- function(
  predictor_variable,
  data_training,
  num_species,
  method) {
  # Determine single vs. multi variable output
  is_single_var = !(predictor_variable %in% c("_abundance", "_composition"))

  if (method == "naive") {
    # Naive fitting doesn't require model
    return(NULL)
  }
  else if (method == "rf") {
    # Random forest prediction
    if (is_single_var) {
      return(fit_rf_singlevar(predictor_variable, data_training, num_species))
    }
    else {
      return(fit_rf_multivar(predictor_variable, data_training, num_species))
    }
  }
  else if (method == "rf_sequential") {
    # Sequential fitting random forest prediction
    stop("Sequential fitting not implemented yet")
  }
  else {
    stop("Invalid fitting method")
  }
}

predict_rf_singlevar <- function(
  predictor_variable,
  rf_model,
  predict_rows) {
  # Predict all values
  values_predicted = predict(
    rf_model, 
    assemblages[predict_rows,]
  )$predictions
  
  return(values_predicted)
}

predict_rf_multivar <- function(
  predictor_variable,
  rf_model,
  predict_rows) {
  # Predict all values
  values_predicted_raw = predict(
    object = rf_model, 
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
    stop("Error, predictor variable is not valid")
  }

  return(values_predicted)
}

predict_naive_singlevar <- function(
  predictor_variable,
  assemblages,
  predict_rows) {
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
  predict_rows) {
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
    stop("Error, predictor variable is not valid")
  }

  return(values_predicted)
}

predict_model <- function(
  predictor_variable,
  rf_model,
  assemblages,
  predict_rows,
  method) {
  # Determine single vs. multi variable output
  is_single_var = !(predictor_variable %in% c("_abundance", "_composition"))

  if (method == "naive") {
    # Naive prediction
    if (is_single_var) {
      return(predict_naive_singlevar(predictor_variable, assemblages, predict_rows))
    }
    else {
      return(predict_naive_multivar(predictor_variable, assemblages, predict_rows))
    }
  }
  else if (method == "rf") {
    # Random forest prediction
    if (is_single_var) {
      return(predict_rf_singlevar(predictor_variable, rf_model, predict_rows))
    }
    else {
      return(predict_rf_multivar(predictor_variable, rf_model, predict_rows))
    }
  }
  else if (method == "rf_sequential") {
    # Sequential fitting random forest prediction
    stop("Sequential predicting not implemented yet")
  }
  else {
    stop("Invalid predicting method")
  }
}

get_ground_truth_values <- function(
  predictor_columns,
  assemblages,
  predict_rows) {
  # Process the ground truth values
  if (predictor_variable == "_composition") {
    values_ground_truth = (assemblages[predict_rows, 1:num_species] > 0)
  }
  else if (predictor_variable %in% c("_abundance", "richness")) {
    values_ground_truth = assemblages[predict_rows, predictor_columns]
  }
  else {
    stop("Error, predictor variable is not valid")
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
    obs = values_observed)
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
    obs = values_observed)
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
  if (is_single_var & 
    (is.factor(assemblages[,predictor_variable]) | 
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
