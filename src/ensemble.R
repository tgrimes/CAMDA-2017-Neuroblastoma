###############################################################################
#
# Bagging
#
###############################################################################
ensemble <- function(results_file, n_bags = 20) {
 
  #results_file <- "output/results/all_results_EFS_2017_07_13_03_02_01"
  load(results_file)
  best_results_data <- all_results$transcripts_introns$model$elnet_model$data
  
  # Remove the transcripts + intron model (stored as all_results[5]).
  temp <- attributes(all_results)
  all_results <- all_results[c(1, 2, 3, 4)]
  for(i in 2:length(temp)) {
    attr(all_results, names(temp)[i]) <- temp[[i]]
  }
  attributes(all_results)
  
  # Fit the ensemble. 
  # Censored observations will be replaced with predicted times from 
  # transcripts + introns model using elastic net (best_results_data).
  ensemble <- fit_ensemble(n_bags, 
                           all_results = all_results, 
                           best_results_data = best_results_data)
  
  # Obtain the predicted values from ensemble on testing set.
  training_index <- ensemble$training_index
  y_hat <- exp(predict.ensemble(ensemble, x = ensemble$x[-training_index, ]))
  
  # Reload results and concatinated it with the ensemble.
  load(results_file)
  y <- all_results[[1]]$pred[[1]]$data$true
  observed <- all_results[[1]]$pred[[1]]$data$observed
  
  all_results$ensemble <- 
    list(pred = list(ensemble = list(data = data.frame(true = y, 
                                                       estimate = y_hat, 
                                                       observed = observed, 
                                                       residual = y_hat - y))),
         ensemble)
  
  #ensemble_coef_summary(ensemble)
  
  survival_type <- attr(all_results, "survival_type")
  save(all_results, file = paste("output/results/with_ensemble", survival_type,
                                 timestamp(), sep = "_"))
}

fit_ensemble <- function(n_bags, 
                         all_results = NULL,
                         results_file = NULL,
                         data_dir = "data/processed/",
                         ensemble_dir = "output/ensemble/",
                         log_dir = "logs/ensemble/",
                         best_results_data = NULL,
                         estimate_censored_times = TRUE) {
  if(is.null(all_results) & is.null(results_file)) {
    stop("In bagging(): Either all_results object or results_file string 
         must be provided.")
    return(NA)
  }
  
  log_file <- paste(log_dir, timestamp(), ".txt", sep = "")
  open_log(log_file)
  cat(format(Sys.time()), "Running `fit_ensemble()`...\n")
  
  seed <- 100 
  set.seed(seed)
  cat("Setting seed to", seed, "\n")

  # Load previous fit:
  if(!is.null(all_results)) {
    cat("Modelling results provided by `all_results`.\n")
  } else if(!is.null(results_file)) {
    cat("Loading results from file", results_file, "\n")
    load(results_file)
  } else {
    cat("Either `all_results` or `results_file` must be provided.\n")
    sink()
    stop("Either `all_results` or `results_file` must be provided.")
  }
  
  print(attributes(all_results))
  
  survival_type <- attr(all_results, "survival_type")
  training_index <- attr(all_results, "training_index")

  # Subset data by high risk for KM assessment.
  load("data/patient_data")
  hr = (patient_data$high_risk == "Y")[training_index] # Boolean vector.
  rm(patient_data)
  
  # Extract observations from results of training. Replace unobserved times with
  # an estimate. Set all observed to true.
  y <- all_results[[1]]$model[[1]]$data$true
  if(estimate_censored_times) {
    y[!best_results_data$observed] <- 
      best_results_data$estimate[!best_results_data$observed]
    observed <- rep(TRUE, length(y))
  }
  
  # Load all predictors used in the various models.
  # The model must spcify x_file_list, which will be used to obtain and index x.
  x_list <- load_x(all_results)
  
  # X matrix containing all predictors.
  x <- x_list$x
  p <- ncol(x)
  
  # Named list (file names) of column indicies of X corresponding to that file:
  x_index <- x_list$x_index
    
  best_model <- rep("", n_bags)
  bag <- vector("list", n_bags)
  cat("\n\n")
  for(i in 1:n_bags) {
    cat(" =================================================================\n")
    cat("Bag", i, "of", n_bags, "\n")
    # 1. Resampling the rows of the data (size n with replacement)
    temp <- resampling(x[training_index, ])
    bag_index <- temp$index
    oob_index <- temp$oob
    cat("oob_index (n = ", 
        length(oob_index), " of ", length(training_index),
        ")", sep = "")
    cat(oob_index, "\n")
    
    pred_results <- NULL
    # Loop over each model result.
    n_results <- length(all_results)
    for(j in 1:n_results) {
      model_list <- get_model_list(all_results[[j]]) # List of functions.
      model_names <- names(model_list)
      model_x_file_list <- all_results[[j]]$x_file_list
      model_x_index <- unlist(x_index[which(names(x_index) %in% model_x_file_list)])
      
      x_name <- gsub(data_dir, "", model_x_file_list)
      x_name <- gsub("_.*", "", x_name)
      x_name <- paste(x_name, collapse = "_")
      
      cat("----------------------------------------------------------------\n",
          "Predictor set ", j, " of ", n_results, " (Bag ", i, "/", n_bags, ")",
          " - ", x_name, ":\n",
          "\t x source     ", paste(model_x_file_list, 
                                    collapse = "\n\t              "), "\n",
          "\t Model list   ", paste(model_names, 
                                    collapse = "\n\t              "), "\n",
          sep = "")
      
      cat("\n")
      results <- run_using(x[training_index, model_x_index], 
                           y, observed, 
                           training_index = bag_index,
                           model_list = model_list,
                           testing_index = oob_index)
      
      # Reset coefficients to align with original X.
      for(k in 1:length(results$pred)) {
        results$pred[[k]]$fit$model_x_index <- model_x_index # Original indexing.
        
        beta_hat <- rep(0, p)
        beta_hat[model_x_index] <- results$pred[[k]]$fit$beta_hat
        results$pred[[k]]$fit$beta_hat <- beta_hat
        
        mean_x <- rep(0, p)
        mean_x[model_x_index] <- results$pred[[k]]$fit$mean_x
        results$pred[[k]]$fit$mean_x <- mean_x
        
        sd_x <- rep(1, p)
        sd_x[model_x_index] <- results$pred[[k]]$fit$sd_x
        results$pred[[k]]$fit$sd_x <- sd_x
      }
      
      names(results$pred) <- paste(x_name, "_", names(results$pred), sep = "")
      
      pred_results <- c(pred_results, results$pred)
      cat("\n")
    } # End for j in 1:n_results
    
    
    # 3. Find best performing model.
    best_model_results <- get_winner(pred_results, hr)
    best_model[i] <- best_model_results$winner
    
    # Print out performance measures and indicate the winner.
    print(best_model_results$ranks)
    cat(best_model[i], "\n")
    
    bag[[i]] <- pred_results[[best_model[i]]]$fit
  } # End for i in 1:n_bags
  
  cat(" =================================================================\n\n")
  cat("Bagging completed\n",
      "Bags:", best_model, "\n")
  
  names(bag) <- best_model
  class(bag) <- "ensemble"
  
  attr(bag, "results_file") <- results_file
  attr(bag, "survival_type") <- survival_type
  attr(bag, "seed") <- seed
  
  ensemble_file_name = paste("ensemble", survival_type, timestamp(), sep = "_")
  ensemble_file <- paste(ensemble_dir, ensemble_file_name, sep = "")
  
  ensemble <- list(bag = bag, 
                   file_dir = ensemble_file, 
                   file_name = ensemble_file_name,
                   x = x,
                   training_index = training_index)
  
  save(ensemble, file = ensemble_file)
  cat("Ensemble results saved to", ensemble_file, "\n")
  
  close_log()
  
  return(ensemble)
}

predict.ensemble <- function(ensemble, x) {
  fit <- ensemble$bag
  n_bags <- length(fit)
  
  y <- rep(0, nrow(x))
  for(i in 1:n_bags) {
    y <- y + predict(fit[[i]], x)
  }
  y <- y / n_bags
  
  return(y)
}

resampling <- function(x, MAX_ITER_OOB_SAMPLING = 100) {
  if(is.null(x) & is.null(n)) stop("Either a matrix x or a number of rows n
                                   must be provided.")
  
  n <- nrow(x)
  epsilon <- 1e-12 # Threshold for variability (detect if var(x[, i]) < epsilon).
  
  # First check that the original data contain variability in columns.
  if(any(apply(x, 2, function(z) var(z) < epsilon))) {
    stop("In resampling: no variability in columns of input data x.")
  }
  
  index <- sample(1:n, n, replace = TRUE)
  oob <- setdiff(1:n, index)
  
  i <- 1
  
  # Flags to check for:
  # 1. Check that the out of bag (OOB) sample is not empty.
  # 2. In OOB sample, columns of x have variability.
  # 3. In index sample, columns of x have variability.
  oob_empty <- (length(oob) == 0)
  no_variance_in_oob <- any(apply(x[oob, ], 2, function(z) var(z) < epsilon))
  no_variance_in_index <- any(apply(x[index, ], 2, function(z) var(z) < epsilon))
  
  while(oob_empty | no_variance_in_oob | no_variance_in_index) {
    # Loop until all conditions are met. Stop after maximum number of iterations.
    i <- i + 1
    if(i > MAX_ITER_OOB_SAMPLING) {
      message <- paste("Resample of data failed. Stopped after", i, "iterations.",
                       "Flags:", oob_empty, no_variance_in_oob, no_variance_in_index)
      stop(message)
    }
    
    index <- sample(1:n, n, replace = TRUE)
    oob <- setdiff(1:n, index)
    
    oob_empty <- length(oob) == 0
    if(!is.null(x)) {
      no_variance_in_oob <- any(apply(x[oob, ], 2, function(z) length(unique(z)) == 1))
      no_variance_in_index <- any(apply(x[index, ], 2, function(z) length(unique(z)) == 1))
    }
  }
  cat("Resampling finished in", i, "iterations.\n")
  return(list(index = index, oob = oob))
}

#' @param results a results object from the output of run_using().
#' @param estimates a p by M matrix of parameter estimates from M different models.
#' @notes Either `estimates` or `results` must be provided.
prune_x <- function(results = NULL, estimates = NULL) {
  # If results are provided, extract the estimates.
  if(!is.null(results)) {
    M <- length(results$model)
    p <- length(results$model[[1]]$fit$beta)
    estimates <- matrix(0, nrow = p, ncol = M)
    for(i in 1:M) {
      estimates[, i] <- results$model[[i]]$fit$beta_hat
    }
  }
  
  # Check that the signs are equal across each model. Zeroes are ignored.
  signs_equal <- apply(estimates, 1, function(x) {
      !(any(sign(x) == 1) & any(sign(x) == -1))
    })
  
  # Check that the estimate is nonzero for at least M - 1 models.
  all_nonzero <- apply(estimates, 1, function(x) {
      sum(abs(x) > 1e-12) >= M - 2
    })
  
  # Check that the estimate is in top q-percentile for each model.
  q <- 0.5
  top_percentile <- TRUE
  for(i in 1:M) {
    top_percentile <- top_percentile & 
      (abs(estimates[, i]) >= quantile(abs(estimates[, i]), q))
  }
  
  # Keep genes that (1) have the same sign across all models, (2) are nonzero
  # for each model excluding the Lasso, (3) are in the top q-percentile for
  # each model.
  keep <- top_percentile
  
  return(keep)
}

get_winner <- function(results, hr) {
  model_names <- names(results)
  n_models <- length(model_names)
  temp <- loss.all(results[[1]]$data, hr)
  ranks <- matrix(0, 
                  nrow = nrow(temp), 
                  ncol = n_models)
  ranks[, 1] <- as.numeric(temp[, 1])
  for(i in 2:n_models) {
    ranks[, i] <- as.numeric(loss.all(results[[i]]$data, hr)[, 1])
  }
  rownames(ranks) <- loss.all(results[[1]]$data, hr)[, 2]
  colnames(ranks) <- model_names
  
  ranking <- matrix(model_names[t(apply(ranks, 1, order))], ncol = n_models)
  all_tied <- which(apply(ranks, 1, function(x) length(unique(x)) == 1))
  if(length(all_tied) > 0) {
    ranking <- ranking[-all_tied, ]
    ranks <- ranks[-all_tied, ]
  }
  
  winner <- ranking[1, 1]
  if(nrow(ranking) > 1) {
    winner <- RankAggreg(ranking, n_models, verbose = FALSE)$top.list[1]
  }
  
  return(list(winner = winner, ranking = ranking, ranks = ranks))
}


load_x <- function(all_results) {
  files <- lapply(all_results, function(result) result$x_file_list)
  files <- unique(unlist(files))
  x <- NULL
  x_index <- vector("list", length(files))
  names(x_index) <- files
  
  for(i in 1:length(files)) {
    temp <- load_RData_file(files[[i]])
    if(is.null(x)) {
      x <- temp
    } else {
      x <- cbind(x, temp)
    }
    x_index[[i]] <- (ncol(x) - ncol(temp) + 1):ncol(x)
  }
  attr(x, "index") <- x_index
  
  return(list(x = x, x_index = x_index))
}

get_model_list <- function(results) {
  model_names <- names(results$model)
  n_models <- length(model_names)
  model_list <- vector("list", n_models)
  names(model_list) <- model_names
  for(i in 1:n_models) {
    model_list[[i]] <- getFunction(model_names[i])
  }
  
  return(model_list)
}

ensemble_coef_summary <- function(ensemble) {
  n <- length(ensemble)
  p <- length(ensemble[[1]]$beta_hat)
  beta_hat <- matrix(0, nrow = n, ncol = p)
  for(i in 1:n) {
    beta_hat[i, ] <- ensemble[[i]]$beta_hat
  }
  colnames(beta_hat) <- colnames(ensemble$x)
  
  all_zeroes <- apply(beta_hat, 2, function(x) all(x == 0))
  x_names[!all_zeroes]
  sign_flips <- apply(beta_hat, 2, function(x) length(unique(sign(x[x != 0]))) > 1)
  beta_hat[, which(sign_flips)[1:min(10, sum(sign_flips))]]
  
  # Index beta hat by order of magnitude. Estimates of zero are ommited.
  index <- order(abs(apply(beta_hat, 2, mean)), decreasing = TRUE)[1:(p - sum(all_zeroes))]
  beta_hat[, index[1:min(10, length(index))]]

  effect_size <- scale(x[training_index, index]) %*% diag(apply(beta_hat[, index], 2, mean))
  effect_size <- apply(abs(effect_size), 2, median)
  by_effect <- order(exp(effect_size), decreasing = TRUE)
  plot(1:length(index), exp(effect_size)[by_effect], type = "l")
  abline(v = 25, col = "gray")
  
  return(data.frame(x = x_names[by_effect], effect_size = effect_size[by_effect]))
}

# # Obtain the various sets of predictors x for each model.
# n_x_sets <- length(all_results)
# x_list <- vector("list", n_x_sets) 
# for(i in 1:n_x_sets) {
#   results <- all_results[[i]]
#   results_data <- results$model$elnet_model$data
#   
#   x_list[[i]]$x <- get_predictors(results$x_file_list)
#   predictors <- which(prune_x(results))
#   #colnames(x[[i]])[predictors]
#   x_list[[i]]$x <- x_list[[i]]$x[training_index, predictors]
#   
#   x_list[[i]]$x_file_list <- results$x_file_list
#   x_list[[i]]$name <- names(all_results)[i]
# }