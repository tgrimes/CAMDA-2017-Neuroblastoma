# library(spls) 

spls_model <- function(y, x, observed = NULL, K = 1:10, eta = seq(0.1, 0.9, 0.1), 
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
  
  # If samples are uncensored, proceed with spls as usual.
  if(is.null(observed) || all(observed == TRUE)) {
    cat("No censoring in the data; SPLS is computed as usual.\n")
    fit_method = "UNCENSORED"
    fit <- fit.spls_model(y, x, K, eta)
    
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
      fit <- impute_data(fit.spls_model, y, x, observed, K = K, eta = eta, ...)
    } else {
      stop('Invalid fit_method for spls model.')
    }
  }
  
  # Fit from model contains beta coefficients.
  # These are concatenated with the standardization values.
  fit <- c(fit, list(mean_y = mean_y, mean_x = mean_x, sd_x = sd_x,
                     fit_method = fit_method))
  class(fit) <- "spls_model"
  return(fit)
}

#' @param y the response variable, assumed to be centered.
#' @param x the explanatory variables, assumed to be standardized.
#' @param v the number of SPLS components to use.
fit.spls_model <- function(y, x, K = NULL, eta = NULL, ...) {
  if(is.null(K)) {
    K <- 1:10 # Default lattice to search on for K if none is provided.
  }
  if(is.null(eta)) {
    eta <- seq(0.1, 0.9, 0.1) # Default lattice to search on for eta.
  }
  
  # Determine K and eta by k-fold cross-validation.
  # cv.spls uses k-fold = 10 by default.
  if(length(eta) > 1 | length(K) > 1) {
    crossval <- cv.spls(x = x, y = y, K = K, eta = eta, plot.it = FALSE, 
                        scale.x = FALSE, scale.y = FALSE) 
    K <- crossval$K.opt
    eta <- crossval$eta.opt
  }
  
  fit <- spls(x = x, y = y, K = K, eta = eta, scale.x = FALSE, scale.y = FALSE)
  
  beta_hat <- coef(fit)
  
  return(list(mu_hat = 0, beta_hat = beta_hat, K = K, eta = eta))
}

#' @param fit fitted model of class 'spls_model'.
#' @param x the explanatory variables. Columns are assumed to have zero mean.
predict.spls_model <- function(fit, x) {
  if(!is.matrix(x)) x <- as.matrix(x)
  
  x <- scale(x, fit$mean_x, fit$sd_x)
  y_hat <- fit$mean_y + fit$mu_hat + x %*% fit$beta_hat
  
  return(y_hat)
}

test.spls_model <- function() {
  n <- 50
  p <- 100
  K = 1
  eta = 0.5
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  true_coef <- rnorm(p) * rbinom(p, 1, 0.2)
  noise <- rnorm(n, 0, 0.2)
  y <- x %*% true_coef + noise
  observed <- rbinom(n, 1, 0.9)
  
  x <- scale(x)
  y <- scale(y)
  
  fit <- spls_model(y, x, observed, K = 1:10)
  cbind(true_coef, fit$beta_hat)
  
  fit <- spls(x, y, K, eta)
  cbind(true_coef, coef(fit))
  
  par(mfrow = c(1, 2))
  newx <- scale(x, fit$meanx, fit$normx)
  pred_y <- newx %*% fit$betahat + mean(y)
  plot(pred_y, pred_y - y, main = mean(y), ylim = c(-3, 3))
  abline(h = 0)
  
  pred_y2 <- predict(fit)
  plot(pred_y2, pred_y2 - y, main = fit$mu, ylim = c(-3, 3))
  abline(h = 0)
}