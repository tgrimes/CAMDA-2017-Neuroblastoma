# library(glmnet) # For glmnet(), which fits lasso.
# library(foreach)
# library(doParallel)

elnet_model <- function(y, x, observed = NULL, alpha = seq(0.1, 0.9, 0.1), lambda = NULL, 
                        scale_x = TRUE, fit_method = "IMPUTE", ...) {
  if(!is.matrix(x)) x <- as.matrix(x)
  
  n <- nrow(x)
  p <- ncol(x)
  one <- matrix(1, 1, n)
  mean_y <- drop(one %*% y) / n
  y <- scale(y, mean_y, FALSE)
  mean_x <- drop(one %*% x) / n
  x <- scale(x, mean_x, FALSE)
  if(scale_x) {
    sd_x <- sqrt(drop(one %*% x^2) / (n - 1))
    if(any(sd_x < .Machine$double.eps)) {
      stop("Some columns of x have zero variance.")
    }
    x <- scale(x, FALSE, sd_x)
  } else {
    sd_x <- rep(1, p)
  }
  
  # If samples are uncensored, proceed with elastic net as usual.
  if(is.null(observed) || all(observed == TRUE)) {
    cat("No censoring in the data; elastic net is computed as usual.\n")
    fit_method = "UNCENSORED"
    fit <- fit.elnet_model(y, x, alpha = alpha, lambda = lambda)
    
  } else {
    # Some observations are censored,so use the imputation method.
    
    # Enforce logical type for variable `observed`.
    if(!is.logical(observed)) {
      if (all(observed == 1 | observed == 0)) {
        observed <- observed == 1 # Convert to boolean.
      } else {
        stop("argument 'observed' should be logical or contain only 1's and 0's.")
      }
    }
    
    # Run imputation algorithm.
    if(fit_method == "IMPUTE") {
      fit <- impute_data(fit.elnet_model, y, x, observed, 
                         alpha = alpha, lambda = lambda, ...)
    } else {
      stop('Invalid fit_method for elnet model.')
    }
  }
  
  # Fit from model contains beta coefficients.
  # These are concatenated with the standardization values.
  fit <- c(fit, list(mean_y = mean_y, mean_x = mean_x, sd_x = sd_x, 
                     fit_method = fit_method))
  class(fit) <- "elnet_model"
  return(fit)
}

#' @param y the response variable - typically log survival times.
#' @param x the explanatory variables - assumed to be centered.
#' @param lambda the penalization parameter.
fit.elnet_model <- function(y, x, alpha = seq(0.1, 0.9, 0.1), lambda = NULL, ...) {
  
  # Determine alpha and lambda by k-fold cross-validation.
  n_alpha <- length(alpha)
  if(n_alpha > 1) {
    # Run in parallel if using Linux machine (won't work in Windows).
    if(.Platform$OS.type == "unix") {
      n_cores <- detectCores()
      cat("Running CV for elastic net on", n_cores, "cores.\n")
      registerDoParallel(n_cores)
      cv_list <- foreach(i = 1:n_alpha) %dopar% {
        crossval <- cv.glmnet(x = x, y = y, alpha = alpha[i], standardize = FALSE,
                              penalty.factor = c(rep(0, 3), rep(1, ncol(x) - 3))) 
        min_index <- which(crossval$cvm == min(crossval$cvm))
        return(list(cvm = crossval$cvm[min_index], lambda = crossval$lambda.min))
      }
      stopImplicitCluster()
      
      cvm <- sapply(cv_list, function(x) x$cvm)
      lambda <- sapply(cv_list, function(x) x$lambda)
      
    # Otherwise, run in sequence.
    } else {
      cvm <- NULL
      lambda <- NULL
      for(i in 1:n_alpha) {
        crossval <- cv.glmnet(x = x, y = y, alpha = alpha[i], standardize = FALSE,
                              penalty.factor = c(rep(0, 3), rep(1, ncol(x) - 3))) 
        min_index <- which(crossval$cvm == min(crossval$cvm))
        cvm <- c(cvm, crossval$cvm[min_index]) 
        lambda <- c(lambda, crossval$lambda.min)
      }
    }
    
    # Find minimum cross-validation error and corresponding alpha and lambda.
    min_index <- which(cvm == min(cvm))
    alpha <- alpha[min_index]
    lambda <- lambda[min_index]
  } else {
    crossval <- cv.glmnet(x = x, y = y, alpha = alpha, standardize = FALSE,
                          penalty.factor = c(rep(0, 3), rep(1, ncol(x) - 3))) 
    lambda <- crossval$lambda.min
  }
  cat("elnet: alpha = ", alpha, ", lambda = ", lambda, 
      ", from 10-fold cv.\n", sep = "")
  fit <- glmnet(x = x, y = y, alpha = alpha, lambda = lambda,
                penalty.factor = c(rep(0, 3), rep(1, ncol(x) - 3)))
  
  
  mu_hat <- coef(fit)[1]
  beta_hat <- coef(fit)[-1]
  
  return(list(mu_hat = mu_hat, beta_hat = beta_hat, alpha = alpha, lambda = lambda))
}

#' @param fit fitted model of class 'elnet_model'.
#' @param x the explanatory variables. 
predict.elnet_model <- function(fit, x) {
  if(!is.matrix(x)) x <- as.matrix(x)
  
  x <- scale(x, fit$mean_x, fit$sd_x)
  y_hat <- fit$mean_y + fit$mu_hat + x %*% fit$beta_hat
  
  return(y_hat)
}

