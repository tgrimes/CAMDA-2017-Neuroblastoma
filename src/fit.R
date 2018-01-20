#' @param x the n by p matrix of predictors.
#' @param y the vector of n observed survival times.
#' @param observed the vector of n logical (or 0/1) indicators for observed event.
#' @param training_index the set of indicies (values 1-n) for the training sample.
#' @param model_list a named list containing the models to use. 
#' @param seed seed used prior to each model fitting. Used for reproducibility.
#' @param testing_index by default is set to the rows not used in training.
#' @param seed for reproducibility, the RNG seed is reset to this value prior
#'   to each model fitting procedure. 
run_using <- function(x, y, observed, training_index,
                      model_list = list(pls_model = pls_model, 
                                        spls_model = spls_model, 
                                        lasso_model = lasso_model, 
                                        elnet_model = elnet_model), 
                      seed = NULL, testing_index = NULL, ...) {
  # Training set:
  x_train <- x[training_index, ]
  y_train <- y[training_index]
  observed_train <- observed[training_index]
  
  # If no testing index is given, use remaining rows after training.
  if(is.null(testing_index)) {
    testing_index <- setdiff(1:nrow(x), training_index)
  }
  # Testing set:
  x_test <- x[testing_index, ]
  y_test <- y[testing_index]
  observed_test <- observed[testing_index]

  rm(x, y) # Remove to free up memory.
  gc()
  
  n_models <- length(model_list)
  model <- vector("list", n_models)
  pred <- vector("list", n_models)
  
  # Fit models.
  for(i in 1:n_models) {
    # For reproducibility, reset the RNG seed at the start of fitting procedure.
    if(!is.null(seed)) set.seed(seed)
    
    cat("Fitting model ", i, " of ", n_models, ".\n", sep = "")
    model[[i]] <- tryCatch(
      fit_on(model_list[[i]], y_train, x_train, observed_train, ...),
      error = function(e) { 
        cat("*****\n")
        cat("Error while fitting ", names(model_list)[i], ":\n", sep = "")
        print(e)
        cat("*****\n")
        return(NA)
      })
    if(!is.na(model)[i]) {
      pred[[i]] <- pred_on(model[[i]]$fit, y_test, x_test, observed_test)
    } 
    cat("\n")
  }
  
  # Remove any models and prediction results that were not able to be fit.
  model[which(is.na(model))] <- NULL
  pred[which(is.na(model))] <- NULL

  names(model) <- names(model_list)
  names(pred) <- names(model_list)
  
  return(list(model = model, pred = pred))
}



fit_on <- function(model.fun, y, x, observed, ...) {
  if(!is.matrix(x)) x <- as.matrix(x)
  
  # Fit the model on log(y).
  fit <- model.fun(log(y), x, observed, ...)
  
  # Compute MSE of fit.
  data <- data.frame(true = y,
                     estimate = exp(predict(fit, x = x)),
                     observed = observed)
  data <- cbind(data, residual = (data$estimate - data$true))
  rmse <- sqrt(sum(data$residual^2 * observed) / sum(observed))
  
  return(c(list(fit = fit), list(rmse = rmse, data = data)))
}



#' @param fit the fitted model - assumed to be on log(y).
pred_on <- function(fit, y, x, observed) {
  if(!is.matrix(x)) x <- as.matrix(x)
  
  # Compute MSPE on observed y.
  data <- data.frame(true = y,
                     estimate = exp(predict(fit, x = x)),
                     observed = observed)
  data <- cbind(data, residual = data$estimate - data$true)
  rmse <- sqrt(sum(data$residual^2 * observed) / sum(observed))
  
  return(list(fit = fit, rmse = rmse, data = data))
}