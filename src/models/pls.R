pls_model <- function(y, x, observed = NULL, K = 1:10,
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
    cat("No censoring in the data; PLS is computed as usual.\n")
    fit_method = "UNCENSORED"
    fit <- fit.pls_model(y, x, K)
    
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
      fit <- impute_data(fit.pls_model, y, x, observed, ...)
    } else {
      stop('Invalid fit_method for pls model.')
    }
  }
  
  # Fit from model contains beta coefficients.
  # These are concatenated with the standardization values.
  fit <- c(fit, list(mean_y = mean_y, mean_x = mean_x, sd_x = sd_x,
                     fit_method = fit_method))
  class(fit) <- "pls_model"
  return(fit)
}

#' @param y the response variable, assumed to be centered.
#' @param x the explanatory variables, assumed to be standardized.
#' @param K the number of PLS components to use.
fit.pls_model <- function(y, x, K = NULL, ...) {
  if(is.null(K)) {
    K <- 1:10 # Default lattice to search on for K if none is provided.
  }
  
  # Determine K by k-fold cross-validation.
  # cv.spls uses k-fold = 10 by default.
  if(length(K) > 1) {
    crossval <- cv.spls(x = x, y = y, K = K, eta = 0, plot.it = FALSE, 
                        scale.x = FALSE, scale.y = FALSE) 
    K <- crossval$K.opt
  }
  
  fit <- spls(x = x, y = y, K = K, eta = 0, scale.x = FALSE, scale.y = FALSE)
  
  beta_hat <- coef(fit)
  
  return(list(mu_hat = 0, beta_hat = beta_hat, K = K))
}

#' pls_model <- function(y, x, observed = NULL, v = 1, 
#'                       center_x = TRUE, scale_x = TRUE, 
#'                       fit_method = "IMPUTE", ...) {
#'   if(!is.matrix(x)) x <- as.matrix(x)
#'   
#'   n <- nrow(x)
#'   p <- ncol(x)
#'   one <- matrix(1, 1, n)
#'   mean_y <- drop(one %*% y) / n
#'   y <- scale(y, mean_y, FALSE)
#'   mean_x <- drop(one %*% x) / n
#'   x <- scale(x, mean_x, FALSE)
#'   if(scale_x) {
#'     sd_x <- sqrt(drop(one %*% x^2) / (n - 1))
#'     if(any(sd_x < .Machine$double.eps)) {
#'       stop("Some columns of x have zero variance.")
#'     }
#'     x <- scale(x, FALSE, sd_x)
#'   } else {
#'     sd_x <- rep(1, p)
#'   }
#'   
#'   # If samples are uncensored, proceed with pls as usual.
#'   if(is.null(observed) || all(observed == TRUE)) {
#'     cat("No censoring in the data; PLS is computed as usual.\n")
#'     fit_method = "UNCENSORED"
#'     fit <- fit.pls_model(y, x, v)
#'     
#'   } else {
#'     # Some observations are censored,so use the imputation method.
#'     
#'     # Enforce logical type for variable `observed`.
#'     if(!is.logical(observed)) {
#'       if (all(observed == 1 | observed == 0)) {
#'         observed <- observed == 1 # Convert to boolean.
#'       } else {
#'         stop("argument 'observed' should be logical or contain only 1's and 0's.")
#'       }
#'     }
#'     
#'     # Run imputation algorithm.
#'     if(fit_method == "IMPUTE") {
#'       fit <- impute_data(fit.pls_model, y, x, observed, v = v, ...)
#'     } else {
#'       stop('Invalid fit_method for pls model.')
#'     }
#'   }
#'   
#'   # Fit from model contains beta coefficients.
#'   # These are concatenated with the standardization values.
#'   fit <- c(fit, list(mean_y = mean_y, mean_x = mean_x, sd_x = sd_x,
#'                      fit_method = fit_method))
#'   class(fit) <- "pls_model"
#'   return(fit)
#' }
#' 
#' #' @param y the response variable - typically log survival times.
#' #' @param x the explanatory variables - assumed to be centered.
#' #' @param v the number of PLS components to use.
#' fit.pls_model <- function(y, x, v = 1, ...) {
#' 
#'   n <- nrow(x)
#'   p <- ncol(x)
#'   
#'   # The PLS coefficients and components.
#'   pls_coef <- array(0, dim = c(p, v)) 
#'   pls_comp <- array(0, dim = c(n, v))
#'   
#'   # The predictors x are assumed to already be standardized (see `pls_model()`).
#'   # The response y is assumed to be centered (required for the PLS computation).
#' 
#'   # Compute the PLS coefficients and components.
#'   for(i in 1:v) {
#'     # The response variable is residual of response after centering.
#'     pls_coef[, i] <- crossprod(x, y) 
#'     pls_coef[, i] <- pls_coef[, i]  / sqrt(sum(pls_coef[, i]^2))
#'     pls_comp[, i] <- x %*% pls_coef[, i]
#'     
#'     x <- tryCatch(x - pls_comp[, i] %*% solve(crossprod(pls_comp[, i])) %*% 
#'                    crossprod(pls_comp[, i], x), 
#'                   error = function(e) { warning(e); return(NA) })
#'     if(is.na(x) && i < v) {
#'       # Only stop if fewer than v components have been computed.
#'       stop(paste("PLS algorithm halted. Only v =", i, 
#'                  "of", v, "components computed."))
#'     }
#'   }
#'   
#'   # Fit a linear model of the original response y using the PLS components.
#'   fit <- lm(y ~ 1 + pls_comp)
#'   mu_hat <- fit$coefficients[1]
#'   comp_beta <- fit$coefficients[-1]
#'   
#'   return(list(mu_hat = mu_hat, beta_hat = pls_coef %*% comp_beta))
#' }

#' @param fit fitted model of class 'pls_model'.
#' @param x the explanatory variables.
predict.pls_model <- function(fit, x) {
  if(!is.matrix(x)) x <- as.matrix(x)
  
  x <- scale(x, fit$mean_x, fit$sd_x)
  y_hat <- fit$mean_y + fit$mu_hat + x %*% fit$beta_hat
  
  return(y_hat)
}

