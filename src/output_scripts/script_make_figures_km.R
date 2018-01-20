#library(survival)
#library(survminer)

plot.km <- function(all_results, cutoff_years = 2, file_name = NULL) {
  if(length(all_results) != 1) stop("Expected `all_results` to have length 1.")
  
  outcome <- attr(all_results, "survival_type")
  training_index <- attr(all_results, "training_index") # Needed for subsetting.
  rate <- 1 / 365 # Multiplier for days. 1/365 gives survival in year.
  rate_label <- "(years)"
  
  # Subset data by high risk, favorable, INSS satage, etc. ----------------------
  load("data/patient_data")
  hr <- patient_data$high_risk[-training_index] == "Y"
  rm(patient_data)
  
  result_name <- names(all_results)
  results <- all_results[[1]]
  method_names <- gsub("_.*", "", names(results$pred))
  n_methods <- length(results$pred)
  
  # Prepare to save plots if file_name is provided.
  if(!is.null(file_name)) {
    # Only 1 method for clinical only model and ensemble.
    if(length(results$pred) == 1) {
      png(file_name, 933, 600, res = 215)
      par(mfcol = c(1, 1))
      
    # Other models should have 4 methods (pls, spls, lasso, and elnet).
    } else {
      png(file_name, 1400, 3600, res = 450)
      par(mfcol = c(n_methods, 1))
    }
  }
  par(mar = c(3, 3, 1.5, 1))
  par(mgp = c(2, 1, 0))
  
  for(i in 1:n_methods) {
    x <- results$pred[[i]]$data$estimate * rate
    y_test <- results$pred[[i]]$data$true * rate
    observed_test <- results$pred[[i]]$data$observed
    
    score <- x > cutoff_years
    sdiff <- tryCatch(survdiff(Surv(y_test[hr], observed_test[hr]) ~ score[hr]),
                      error = function(e) NULL)
    if(is.null(sdiff)) next
    pvalue <- 1 - pchisq(sdiff$chisq, 1)
    #if(pvalue > 0.059 & k > 1) next
    
    title <- paste("LPS", cutoff_years, sep = "")
    y_lab <- ifelse(outcome == "OS", "Overall Survival", "Event-Free Survival")
    if(result_name == method_names[i]) {
      method_name <- method_names[i]
    } else {
      method_name <- paste(result_name, method_names[i], sep = "_")
    }
    
    pred_survival = c("< 2 years", ">= 2 years")[1 + (score * 1)]
    df <- data.frame(y = y_test[hr], 
                     observed = observed_test[hr],
                     pred_survival = pred_survival[hr])
    sfit <- survfit(Surv(y, observed) ~ pred_survival, data = df)
    
    # main <- paste("HR and ", title, sep = "")
    
    xname <- result_name
    if(length(grep("clinical", xname, ignore.case = TRUE)) == 1) {
      main <- "clinical"
    } else if(length(grep("ensemble", xname, ignore.case = TRUE)) == 1) {
      main <- "ensemble"
    } else {
      xname <- gsub("_", " ", xname)
      main <- paste(xname, method_names[i], sep = " - ")
    }
    
    plot(sfit, main = main, 
         #conf.int = TRUE, 
         mark.time = TRUE,
         xlab = " Time (in years)", ylab = y_lab,
         col = c("orange", "blue"))
    text(x = 1.5, y = 0.05, labels = paste("p =", round(pvalue, 5)))
  }
  
  # Conclude saving plots if file_name was provided.
  if(!is.null(file_name)) {
    dev.off()
  }
}

make_km_plots <- function(results_list, figures_dir, cutoff_years = 2) {
  for(results_file in results_list) {
    load(results_file)
    outcome <- attr(all_results, "survival_type") 
    
    # File name for a given measure.
    get_file_name <- function(method_name) {
      paste(figures_dir, "km/", 
            outcome, "_LPS_", cutoff_years, "_", method_name, ".png", sep = "")
    }
    
    # Create figure for each model.
    for(i in 1:length(all_results)) {
      temp <- all_results
      temp[-i] <- NULL # Remove all but the ith result.
      plot.km(temp, cutoff_years, get_file_name(names(temp)))
    }
  }
}



