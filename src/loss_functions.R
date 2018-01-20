######
#
# Loss functions
#
######

loss.rmse <- function(df, rate = 1 / 365) {
  name <- "RMSE"
  fn <- function(x) {
    km <- survfit(Surv(x$resid, 1 - x$observed) ~ 1, type = "kaplan-meier")
    sorted_resid <- km$time
    sorted_observed <- km$n.censor
    weights <- c(1, km$surv[1:(length(km$surv) - 1)]) # S^c(T^c-)
    sqrt(mean(sorted_observed * sorted_resid^2 / weights))
  }
  x <- data.frame(resid = df$residual * rate,
                  observed = df$observed)
  return(list(loss = fn(x), name = name))
}

loss.lps5 <- function(df, hr, rate = 1 / 365) {
  name <- "LPS5"
  fn <- function(x) {
    test <- survdiff(Surv(x$true, x$observed) ~ x$score)
    log_pvalue <- log(1 - pchisq(test$chisq, 1))
    # Return the test statistic instead; log_pvalues are often -Inf and uncomparable.
    return(-test$chisq)
  }
  x <- data.frame(score = df$estimate[hr] * rate > 5,
                  true = df$true[hr] * rate,
                  observed = df$observed[hr])
  # If survdiff fails, return -chisq = 0 for the test statistic.
  loss <- tryCatch(fn(x), error = function(e) { print(e); 0 })
  return(list(loss = loss, name = name))
}

loss.lps2 <- function(df, hr, rate = 1 / 365) {
  name <- "LPS2"
  fn <- function(x) {
    test <- survdiff(Surv(x$true, x$observed) ~ x$score)
    log_pvalue <- log(1 - pchisq(test$chisq, 1))
    # Return the test statistic instead; log_pvalues are often -Inf and uncomparable.
    return(-test$chisq)
  }
  x <- data.frame(score = df$estimate[hr] * rate > 2,
                  true = df$true[hr] * rate,
                  observed = df$observed[hr])
  # If survdiff fails, return -chisq = 0 for the test statistic.
  loss <- tryCatch(fn(x), error = function(e) { print(e); 0 })
  return(list(loss = loss, name = name))
}

loss.all <- function(df, hr = rep(T, length(df$estimate)), rate = 1 / 365) {
  measures <- rbind()
  return(rbind(loss.rmse(df, rate),
               loss.lps5(df, hr, rate),
               loss.lps2(df, hr, rate)))
}

loss.all2 <- function(y_hat, y, observed, 
                      hr = rep(T, length(y)), rate = 1 / 365) {
  df <- data.frame(true = y,
                     estimate = y_hat,
                     observed = observed,
                     residual = y_hat - y)
  return(loss.all(df, rate))
}
