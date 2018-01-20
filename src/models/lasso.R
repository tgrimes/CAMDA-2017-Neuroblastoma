# library(glmnet) # For glmnet(), which fits lasso.
# library(foreach)
# library(doParallel)

lasso_model <- function(y, x, observed = NULL, lambda = NULL, 
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
  
  # If samples are uncensored, proceed with lasso as usual.
  if(is.null(observed) || all(observed == TRUE)) {
    cat("No censoring in the data; lasso is computed as usual.\n")
    fit_method = "UNCENSORED"
    fit <- fit.lasso_model(y, x, lambda)
    
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
      fit <- impute_data(fit.lasso_model, y, x, observed, lambda = lambda, ...)
    } else {
      stop('Invalid fit_method for lasso model.')
    }
  }
  
  # Fit from model contains beta coefficients.
  # These are concatenated with the standardization values.
  fit <- c(fit, list(mean_y = mean_y, mean_x = mean_x, sd_x = sd_x,
                     fit_method = fit_method))
  class(fit) <- "lasso_model"
  return(fit)
}

#' @param y the response variable - typically log survival times.
#' @param x the explanatory variables - assumed to be centered.
#' @param lambda the penalization parameter.
fit.lasso_model <- function(y, x, lambda = NULL, ...) {
  if(any(is.infinite(y))) stop("which generated y's are zero:", which(is.infinite(y)))
  
  # # Run in parallel.
  # n_cores <- detectCores()
  # cat("Running lasso cv.glmnet() on", n_cores, "cores.\n")
  # registerDoParallel(n_cores)
  
  # Determine lambda by k-fold cross-validation.
  # cv.glmnet uses k = 10 and alpha = 1 (Lasso) by default.
  if((is.null(lambda) | length(lambda) > 1)) {
    crossval <- cv.glmnet(x = x, y = y,
                          penalty.factor = c(rep(0, 3), rep(1, ncol(x) - 3))) 
    lambda <- crossval$lambda.min
  }
  cat("lasso: lambda = ", lambda, ", from 10-fold cv.\n", sep = "")
  fit <- glmnet(x = x, y = y, alpha = 1, lambda = lambda,
                penalty.factor = c(rep(0, 3), rep(1, ncol(x) - 3)))
  
  
  # stopImplicitCluster()
  
  
  mu_hat <- coef(fit)[1]
  beta_hat <- coef(fit)[-1]
  
  return(list(mu_hat = mu_hat, beta_hat = beta_hat, lambda = lambda))
}

#' @param fit fitted model of class 'lasso_model'.
#' @param x the explanatory variables.
predict.lasso_model <- function(fit, x) {
  if(!is.matrix(x)) x <- as.matrix(x)
  
  x <- scale(x, fit$mean_x, fit$sd_x)
  y_hat <- fit$mean_y + fit$mu_hat + x %*% fit$beta_hat
  
  return(y_hat)
}

