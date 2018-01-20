#' @param y the response variable - typically log survival times.
#' @param observed Boolean vector (or 1's and 0's) of length n indicating if 
#'        failure was observed (TRUE or 1 means no censoring).
#' @param x n by p matrix of predictors
#' @param model the AFT model to be fit as a function of uncensored times and x.
#' @param impute_iter (tuning) the number of iterations to perform.
#' @param impute_datasets (tuning) the number of datasets to be imputed.
#' @param ... other parameters to be passed into model.
impute_data <- function(fit.fun, y, x, observed,
                        impute_iter = 5, impute_datasets = 5, 
                        verbose = TRUE, plots = FALSE,...) {
  par(mfrow = c(2, 3))
  
  # Enforce logical type for variable `observed`.
  if(!is.logical(observed)) {
    if (all(observed == 1 | observed == 0)) {
      observed <- observed == 1 # Convert to boolean.
    } else {
      stop("argument 'observed' should be logical or contain only 1's and 0's.")
    }
  }
  
  # If all are observed, do not impute. Use fit.fun and return now.
  if(all(observed)) {
    warning("All observations are observed. Returning fit.fun(y, x, ...).")
    return(fit.fun(y, x, ...))
  }
  
  # Make sure explanatory variables are in matrix form.
  if(!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  n <- nrow(x) # Number of observations.
  p <- ncol(x) # Number of parameters.
  n_censored <- sum(!observed)
  
  # Initialization: Fit the AFT model to uncensored observations.
  if(verbose) cat("Initialization: Fitting AFT model to uncensored observations.\n")
  fit <- fit.fun(y[observed], x[observed, ], ...)

  mu_hat <- fit$mu_hat
  beta_hat <- fit$beta_hat
  
  print(summary(beta_hat))
  # Using the initial estimate of beta, compute the resdiuals.
  residuals <- y - (mu_hat + x %*% beta_hat)
  if(plots) hist(residuals, main = "Initial")
  
  # Create generator for residuals based on Kaplan-Meier estimate of the 
  # marginal distribution of the model error using the product limit formula.
  rresid <- kaplan_meier_conditional(residuals, observed)$generator

  if(verbose) cat("Begin imputation iterations.\n")
  for(j in 1:impute_iter) {
    if(verbose) cat(" - iteration ", j, " of ", impute_iter, ".\n", sep = "")
    
    # Generate impute_datasets new data sets.
    # Fit the model on each generated data set. Average the estimated betas.
    generated_mu_hats <- matrix(0, 1, impute_datasets)
    generated_beta_hats <- matrix(0, p, impute_datasets)

    for(i in 1:impute_datasets) {
      if(verbose) cat("\t Imputed dataset ", i, " of ", impute_datasets, ".\n", sep = "")
      generated_y <- y
      # resid <- rresid(sum(!observed))
      resid <- rresid(residuals[!observed])
      # hist(resid, main = paste("sampled residuals", i))

      generated_y[!observed] <- mu_hat + x[!observed, ] %*% beta_hat + resid
      #generated_y[!observed] <- y[!observed] + resid
      
      # Fit model based on generated censoring.
      fit <- fit.fun(generated_y, x, ...)
      generated_mu_hats[i] <- fit$mu_hat
      generated_beta_hats[, i] <- fit$beta_hat
    }
    
    mu_hat <- mean(generated_mu_hats)
    beta_hat <- apply(generated_beta_hats, 1, mean)

    # Using the updated estimate of mu and beta, recompute the resdiuals.
    residuals <- y - (mu_hat + x %*% beta_hat)
    if(plots) hist(residuals, main = paste("After iteration", j))
    
    # Using the updated estimate of beta, compute the resdiuals and update
    # the Kaplan-Meier estimate. 
    rresid <- kaplan_meier_conditional(residuals, observed)$generator
    
    # par(mfrow = c(3, 2))
    # hist(residuals, main = "Original residuals")
    
    # TRY PARAMETRIC VERSION (USING SURVIVAL)
    
    # Create generator for residuals based on Kaplan-Meier estimate of the 
    # marginal distribution of the model error using the product limit formula.
    # rresid <- kaplan_meier(residuals, observed)$generator
  }
  
  par(mfrow = c(1, 1))
  return(list(mu_hat = mu_hat, beta_hat = beta_hat))
}

kaplan_meier_conditional <- function(residuals, observed) {
  # library(survival)
  # Force largest observation to be uncensored.
  # observed[order(residuals, decreasing = TRUE)[1]] <- 1
  # NULL model rMSPE = 2342.80
  km <- survfit(Surv(residuals, observed) ~ 1, type = "kaplan-meier")
  F_values <- sort((1 - km$surv) / (1 - min(km$surv)))
  original_residuals <- sort(km$time)
  
  generator <- function(condition_residuals) {
    sampled_residuals <- sapply(condition_residuals, function(residual) {
      # Find the index of this residual.
      # residual = 0
      # index <- which(original_residuals >= residual)
      # index <- index[1]
      index <- which(original_residuals == residual)
      if(length(index) < 1) {
        warning(paste("index was not found. residual =", residual))
        index <- 1
      }
      # Sample new residual conditioned on the current residual.
      lower_F_value <- F_values[index]
      if(is.na(lower_F_value)) {
        cat("index =", index, "residual =", residual, "\n",
            "\t F_values = (", length(F_values), ")", F_values[(length(F_values) - 10):(length(F_values))], "\n",
            "\t residuals = (", length(original_residuals), ")", original_residuals[(length(original_residuals) - 10):length(original_residuals)], "\n")
      }
      u <- runif(1, lower_F_value, 1)
      upper <- which(u <= F_values)
      index <- ifelse(length(upper) > 0, upper[1], km$n)
      return(original_residuals[index])
    })
    
    return(sampled_residuals)
  }
  
  return(list(generator = generator))
}








parametric_version <- function(residuals, observed) {
  logLikCensoredFun <- function(params) {
    mu <- params[1]
    sigma <- params[2]
    
    if(sigma < 0) return(NA)
    
    sum(ifelse(observed, 
               log(dnorm(residuals, mu, sigma)),
               log(1 - pnorm(residuals, mu, sigma))))
  }
  cat(sd(residuals), "\t", sd(residuals[observed]), "\n")
  # library(maxLik)
  fit <- maxLik(logLik = logLikCensoredFun,
                start = c(mean = 0, sd = 1))
  
  mu_hat <- fit$estimate[1]
  sigma_hat <- fit$estimate[2]
  
  cat(mu_hat, "\t", sigma_hat, "\n")
  # library(truncnorm)
  generator <- function(condition_residuals) {
    sampled_residuals <- sapply(condition_residuals, function(residual) {
      rtruncnorm(1, residual, Inf, mu_hat, sigma_hat)
    })
    
    return(sampled_residuals)
  }
  
  return(list(generator = generator, mu_hat, sigma_hat))
}

parametric_version.test <- function(n = 249, censor_rate = 0.6) {
  par(mfrow = c(1, 3))
  
  residuals <- rnorm(n)
  residuals <- sort(residuals)
  delta <- rbinom(n, 1, prob = 1 - censor_rate)
  observed <- delta == 1
  print(Surv(residuals, observed))
  
  # For censored residuals, subtract a random amount from exponential dist.
  residuals <- sapply(1:n, function(i) {
    if(delta[i] == 0) {
      return(residuals[i] - rexp(1, 1))
    } else if(delta[i] == 1) {
      return(residuals[i])
    }
  })
  
  hist(residuals, main = paste("Sampled residuals\n",
                               "mean =", round(mean(residuals), 2)))
  points(x = residuals[!observed], y = rep(0, sum(!observed)), col = "red")
  
  rresid <- parametric_version(residuals, observed)
  
  plot(sort(residuals), rresid$F_values, type = "l", xlim = c(-3, 3), ylim = c(0, 1),
       main = paste("Kaplan_meier()$F_values\n n =", 
                    n, "censor rate =", censor_rate), 
       xlab = "survival time", ylab = "F")
  lines(seq(-3, 3, 0.01), pnorm(seq(-3, 3, 0.01)), col = "red")
  
  resid_sample <- rresid$generator(rnorm(10000))
  hist(resid_sample, freq = FALSE, xlim = c(-3, 3),
       main = paste("kaplan_meier()$generator\n n =", 
                    n, "censor rate =", censor_rate), 
       xlab = "survival time")
  lines(seq(-3, 3, 0.01), dnorm(seq(-3, 3, 0.01)), col = "red")
}

kaplan_meier <- function(residuals, observed) {
  cat("kaplan_meier\n")
  
  # library(survival)
  # Force largest observation to be uncensored.
  # observed[order(residuals, decreasing = TRUE)[1]] <- 1
  # NULL model rMSPE = 2342.80
  km <- survfit(Surv(residuals, observed) ~ 1, type = "kaplan-meier")
  km$surv <- 1 - (1 - km$surv) / (1 - min(km$surv))
  
  generator <- function(n) {
    if(n == 0) {
      warning("Attempted to generate 0 observations.")
      return(NULL)
    }
    u <- runif(n)
    indicies <- sapply(u, function(x) {
      upper <- which(x <= 1 - km$surv)
      val <- ifelse(length(upper) > 0, upper[1], km$n)
      return(val)
    })
    sampled_residuals <- km$time[indicies]
    return(sampled_residuals)
  }
  
  return(list(generator = generator, F_values = 1 - km$surv, km = km))
}

kaplan_meier.test <- function(n = 250, censor_rate = 0.6) {
  par(mfrow = c(1, 3))
  
  residuals <- rnorm(n)
  residuals <- sort(residuals)
  delta <- rbinom(n, 1, prob = 1 - censor_rate)
  observed <- delta == 1
  print(Surv(residuals, observed))
  
  # For censored residuals, subtract a random amount from exponential dist.
  residuals <- sapply(1:n, function(i) {
    if(delta[i] == 0) {
      return(residuals[i] - rexp(1, 1))
    } else if(delta[i] == 1) {
      return(residuals[i])
    }
  })
  
  hist(residuals, main = paste("Sampled residuals\n",
                               "mean =", round(mean(residuals), 2)))
  points(x = residuals[!observed], y = rep(0, sum(!observed)), col = "red")
  
  rresid <- kaplan_meier_conditional(residuals, observed)
  
  plot(sort(residuals), rresid$F_values, type = "l", xlim = c(-3, 3), ylim = c(0, 1),
       main = paste("Kaplan_meier()$F_values\n n =", 
                    n, "censor rate =", censor_rate), 
       xlab = "survival time", ylab = "F")
  lines(seq(-3, 3, 0.01), pnorm(seq(-3, 3, 0.01)), col = "red")
  
  resid_sample <- rresid$generator(sample(residuals, 10000, replace = TRUE))
  print(summary(resid_sample))
  hist(resid_sample, freq = FALSE, xlim = c(-3, 3),
       main = paste("kaplan_meier()$generator\n n =", 
                    n, "censor rate =", censor_rate), 
       xlab = "survival time")
  lines(seq(-3, 3, 0.01), dnorm(seq(-3, 3, 0.01)), col = "red")
}





