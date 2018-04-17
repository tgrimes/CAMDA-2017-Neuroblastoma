#library(survival)
#library(boot)

get_measures <- function(all_results, B = 1000, seed = 12345) {
  set.seed(seed)
  # Load results ----------------------------------------------------------------
  outcome <- attr(all_results, "survival_type")
  training_index <- attr(all_results, "training_index") # Boolean vector.
  rate <- 1 / 365 # Multiplier for days. 1/365 gives survival in year.
  rate_label <- "(years)"
  
  # Subset data by high risk, favorable, INSS satage, etc. ----------------------
  load("data/patient_data")
  hr = (patient_data$high_risk == "Y")[-training_index]
  rm(patient_data)
  
  result_names <- names(all_results)
  n_results <- length(result_names)
  
  # Table to store results.
  measures <- NULL
  
  # Function to estimate performance for given measure. 
  get_measures.inner <- function(name, fn, x) {
    bs <- tryCatch(boot(x, fn, B),
                   error = function(e) NULL)
    if(is.null(bs)) {
      warning("boot failed. Setting confidence limits to NA.")
      ci <- c(NA, NA)
      est <- tryCatch(fn(x), error = function(e) NA)
    } else {
      ci <- tryCatch(boot.ci(bs, type = "perc")$percent[4:5],
                     error = function(e) { warning(e); c(NA, NA) })
      if(is.null(ci)) {
        ci <- c(NA, NA) # boot.ci() may fail with warning and return null.
      }
      est <- bs$t0
    }
    
    return(data.frame(data = result_names[k],
                      model = model_names[i],
                      measure = name, 
                      estimate = est[1], 
                      lower = ci[1], 
                      upper = ci[2], stringsAsFactors = FALSE))
  }
  
  for(k in 1:n_results) {
    results <- all_results[[k]]
    
    model_names <- names(results$pred)
    n_models <- length(model_names)
    
    for(i in 1:n_models) {
      cat("iter k =", k, ", i =", i, "\n")
      
      # RMSE (Inverse Probability of Censoring Weighting IPCS)
      name <- "RMSE"
      fn <- function(x, index) {
        km <- survfit(Surv(x$resid[index], 1 - x$observed[index]) ~ 1, 
                      type = "kaplan-meier")
        sorted_resid <- km$time
        sorted_observed <- km$n.censor
        weights <- c(1, km$surv[1:(length(km$surv) - 1)])
        sqrt(mean(sorted_observed * sorted_resid^2 / weights))
      }
      x <- data.frame(resid = results$pred[[i]]$data$residual * rate,
                      observed = (results$pred[[i]]$data$observed))
      measures <- rbind(measures, get_measures.inner(name, fn, x))
      
      
      # Harrell's C-index:
      name <- "c-index"
      fn <- function(x, index) {
        x <- x[index, ]
        # By default, scores in the opposite direction as true times are concordant.
        # Use 1 - c-index.
        results <- survConcordance(Surv(x$true, x$observed) ~ x$score)
        if(all(index == 1:nrow(x))) {
          
        }
        return(1 - results$concordance)
      }
      x <- data.frame(score = results$pred[[i]]$data$estimate,
                      true = results$pred[[i]]$data$true,
                      observed = (results$pred[[i]]$data$observed))
      measures <- rbind(measures, get_measures.inner(name, fn, x))
      # fit <- survConcordance(Surv(x$true, x$observed) ~ x$score)
      # measures <- rbind(measures,
      #                   data.frame(data = result_names[k],
      #                              model = model_names[i],
      #                              measure = name, 
      #                              estimate = 1 - fit$concordance, 
      #                              lower = 1 - (fit$concordance + 1.96 * fit$std.err), 
      #                              upper = 1 - (fit$concordance - 1.96 * fit$std.err), 
      #                              stringsAsFactors = FALSE))
      
      
      # Classification fit using 5 year cutoff.
      name <- "LPS5"
      fn <- function(x, index) {
        test <- tryCatch(survdiff(Surv(x$true[index], x$observed[index]) ~ 
                                    x$score[index]),
                         error = function(e) NA)
        val <- ifelse(!is.na(test), test$chisq, 0)
        #pvalue <- 1 - pchisq(test$chisq, 1)
        return(val)
      }
      x <- data.frame(score = results$pred[[i]]$data$estimate * rate > 5,
                      true = results$pred[[i]]$data$true * rate,
                      observed = (results$pred[[i]]$data$observed))
      measures <- rbind(measures, get_measures.inner(name, fn, x[hr, ]))
      
      
      # Classification fit using 2 year cutoff.
      name <- "LPS2"
      fn <- function(x, index) {
        test <- tryCatch(survdiff(Surv(x$true[index], x$observed[index]) ~ 
                                    x$score[index]),
                         error = function(e) NA)
        val <- ifelse(!is.na(test), test$chisq, 0)
        #pvalue <- 1 - pchisq(test$chisq, 1)
        return(val)
      }
      x <- data.frame(score = results$pred[[i]]$data$estimate * rate > 2,
                      true = results$pred[[i]]$data$true * rate,
                      observed = (results$pred[[i]]$data$observed))
      measures <- rbind(measures, get_measures.inner(name, fn, x[hr, ]))
    }
  }
  
  return(measures)
}


plot.measures <- function(measures, 
                          name,  
                          file_name = NULL,
                          outcome = c("OS", "EFS")) {
  outcome <- outcome[1]
  index <- which((measures$measure == name))
  if(length(index) == 0) {
    warning(paste("measure", name, "not found."))
    return()
  }
  measures <- measures[index, ]
  n <- length(index)
  if(length(grep("LPS", name)) > 0) {
    measures$lower <- log10(1 - pchisq(measures$lower, 1))
    measures$upper <- log10(1 - pchisq(measures$upper, 1))
    measures$estimate <- log10(1 - pchisq(measures$estimate, 1))
    y_range <- c(min(measures$upper), max(measures$lower))
    #y_range <- c(min(measures$lower), max(measures$upper))
    y_name <- "log10 p-value"
    
  } else {
    y_range <- c(min(0, min(measures$lower)), max(measures$upper, 0))
    y_name <- name
  }
  color <- c("black", rep(c("blue", "red", "orange", "green"), 4), "gray")
  pch <- c(20, rep(c(1, 4, 0, 5), each = 4), 8)
  
  
  # Prepare to save plots if file_name is provided.
  if(!is.null(file_name)) {
    png(file_name,  1600, 600, res = 180)
    par(mfrow = c(1, 1))
    par(mar = c(3, 3, 1, 1))
    par(mgp = c(2, 1, 0))
  }
  
  main <- paste(ifelse(outcome == "OS", "Overall survival", "Event-free survival"),
                name, sep = " - ")
  x_lab <- "Model"
  x_labels <- c("Null", 
                paste("G", 1:4, sep = "-"), 
                paste("T", 1:4, sep = "-"),
                paste("I", 1:4, sep = "-"),
                paste("TI", 1:4, sep = "-"),
                paste("Ens"))
  plot(x = 1:n, y = measures$estimate, pch = 19, col = "black", #pch = pch, col = color, 
       ylim = y_range, ylab = y_name, xlab = "Model", xaxt = "n",
       main = main)
  axis(1, at = 1:length(x_labels), labels = x_labels)
  for(i in 1:n) {
    lines(x = rep(i, 2), y = c(measures$lower[i], measures$upper[i]), 
          lwd = 3, col = "black") #col = color[i])
    points(x = rep(i, 2), y = c(measures$lower[i], measures$upper[i]),
           pch = "-", cex = 2, col = "black") #col = color[i])
  }
  abline(h = 0)
  
  # Conclude saving plots if file_name was provided.
  if(!is.null(file_name)) {
    dev.off()
  }
}

make_performance_plots <- function(results_list, figures_dir, B = 1000) {
  for(results_file in results_list) {
    load(results_file)
    outcome <- attr(all_results, "survival_type")
    print(outcome)
    measures <- get_measures(all_results, B = B)
    
    # File name for a given measure.
    get_file_name <- function(measure) {
      paste(figures_dir, "measures/", 
            outcome, "_", measure, ".png", sep = "")
    }
    
    # Create figure for each measure.
    for(measure in unique(measures$measure)) {
      plot.measures(measures, measure, get_file_name(measure), outcome)
    }
  }
}




