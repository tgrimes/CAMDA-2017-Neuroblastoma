simple_model <- function(y, x, observed = NULL, scale_x = TRUE, 
                         fit_method = "IMPUTE", ...) {
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
  
  # If samples are uncensored, proceed with simple model as usual.
  if(is.null(observed) || all(observed == TRUE)) {
    cat("No censoring in the data; linear regression is computed as usual.\n")
    fit_method = "UNCENSORED"
    fit <- fit.simple_model(y, x)
    
  } else {
    # Some observations are censored, so use the imputation method.
    
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
      fit <- impute_data(fit.simple_model, y, x, observed, ...)
    } else {
      stop('Invalid fit_method for simple model.')
    }
  }
  
  # Fit from model contains beta coefficients.
  # These are concatenated with the standardization values.
  fit <- c(fit, list(mean_y = mean_y, mean_x = mean_x, sd_x = sd_x, 
                     fit_method = fit_method))
  class(fit) <- "simple_model"
  return(fit)
}

#' @param y the response variable - typically log survival times.
#' @param x the explanatory variables - assumed to be centered.
fit.simple_model <- function(y, x) {
  fit <- lm(y ~ 1 + x)
  mu_hat <- fit$coefficients[1]
  beta_hat <- fit$coefficients[-1]
  
  return(list(mu_hat = mu_hat, beta_hat = beta_hat))
}

#' @param fit fitted model of class 'simple_model'.
#' @param x the explanatory variables.
predict.simple_model <- function(fit, x) {
  if(!is.matrix(x)) x <- as.matrix(x)
  
  x <- scale(x, fit$mean_x, fit$sd_x)
  y_hat <- fit$mean_y + fit$mu_hat + x %*% fit$beta_hat
  
  return(y_hat)
}

