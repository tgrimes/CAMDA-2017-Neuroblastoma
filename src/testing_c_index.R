times <- rnorm(100)
x <- data.frame(score = times + runif(100, -1, 1),
                true = times,
                observed = rbinom(100, 1, 0.4))

get_c_index <- function(actual, predicted, delta) {
  n <- length(actual)
  numerator <- 0 # Sum i != j Delta_i*I(pred_i <= pred_j)*I(true_i <= true_j)
  denominator <- 0 # # Sum i != j Delta_i*I(pred_i <= pred_j)
  for(i in 1:n) {
    if(delta[i] == 0) next # Skip if smaller value is censored.
    for(j in 1:n) {
      if(i == j) next # Sum over i != j.
      indicator_1 <- (predicted[i] < predicted[j])
      indicator_2 <- (actual[i] < actual[j])
      numerator <- numerator + indicator_1 * indicator_2
      denominator <- denominator + indicator_1
    }
  }
  return(numerator / denominator)
}
fit1 <- get_c_index(x$score, x$true, x$observed)

library(survcomp)
fit2_all <- concordance.index(x$score, surv.time = x$true, surv.event = x$observed)
fit2 <- 1 - fit2_all$c.index

library(survival)
fit3_all <- survConcordance(Surv(x$true, x$observed) ~ x$score)
fit3 <- 1 - fit3_all$concordance
se3 <- fit3_all$std.err

fit3 - 1.96 * se3
fit3 + 1.96 * se3

fit4_all <- rcorrcens(Surv(x$true, x$observed) ~ x$score)
fit4 <- fit4_all[, "C"]
se4 <- fit4_all[, "SD"]

fit5_all <- rcorr.cens(x$score, Surv(x$true, x$observed))
fit5 <- fit5_all["C Index"]
se5 <- fit5_all["S.D."]

microbenchmark::microbenchmark(
  # get_c_index(x$score, x$true, x$observed),
  concordance.index(x$score, surv.time = x$true, surv.event = x$observed),
  survConcordance(Surv(x$true, x$observed) ~ x$score),
  rcorrcens(Surv(x$true, x$observed) ~ x$score),
  times = 100
)
