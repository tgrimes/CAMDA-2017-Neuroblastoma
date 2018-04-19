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
      main <- "null model"
    } else if(length(grep("ensemble", xname, ignore.case = TRUE)) == 1) {
      main <- "ensemble"
    } else {
      model_name_list <- c("genes", "transcripts", "introns", "transcripts_introns")
      model_symbol_list <- c("G", "T", "I", "TI")
      
      method_name_list <- c("pls", "spls", "lasso", "elnet")
      method_symbol_list <- c("1", "2", "3", "4")
      
      main <- paste(model_symbol_list[which(model_name_list %in% xname)], 
                    method_symbol_list[which(method_name_list %in% method_names[i])], 
                    sep = "-")
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




plot.km.all <- function(all_results, cutoff_years = 2, file_name = NULL,
                        use_hr = FALSE) {
  if(length(all_results) != 1) stop("Expected `all_results` to have length 1.")
  
  outcome <- attr(all_results, "survival_type")
  training_index <- attr(all_results, "training_index") # Needed for subsetting.
  rate <- 1 / 365 # Multiplier for days. 1/365 gives survival in year.
  rate_label <- "(years)"
  
  # Classify patients by high risk, favorable, INSS satage, etc. ----------------------
  load("data/patient_data")
  hr <- patient_data$high_risk[-training_index] == "Y"
  rm(patient_data)
  
  result_name <- names(all_results)
  results <- all_results[[1]]
  method_names <- gsub("_.*", "", names(results$pred))
  n_methods <- length(results$pred)
  
  # Prepare to save plots if file_name is provided.
  if(!is.null(file_name)) {
    png(gsub("/(OS|EFS)[^/]*.png", "/hr_only.png", file_name), 
        1600, 1200, res = 350)
    par(mfcol = c(1, 1))
  }
  par(mar = c(3, 3, 1.5, 1))
  par(mgp = c(2, 1, 0))
  # Plot KM curves for subset only.
  df <- data.frame(y = results$pred[[1]]$data$true * rate, 
                   observed = results$pred[[1]]$data$observed,
                   hr = hr)
  sfit <- survfit(Surv(y, observed) ~ hr, data = df)
  title <- paste0("hr (", sfit$strata[2], ":", sfit$strata[1], ")")
  y_lab <- ifelse(outcome == "OS", "Overall Survival", "Event-Free Survival")
  plot(sfit, main = title, 
       #conf.int = TRUE, 
       mark.time = TRUE,
       xmax = 8,
       xlab = "Time (in years)", 
       ylab = y_lab,
       col = c("black", "red"))
  if(!is.null(file_name)) {
    dev.off()
  }
  
  # Prepare to save plots if file_name is provided.
  if(!is.null(file_name)) {
    if(use_hr) {
      file_name <- gsub("/[^/]*.png", 
                        paste0("/", outcome, "_hr_only.png"), 
                        file_name)
    } else {
      file_name <- file_name
    }
    # Only 1 method for clinical only model and ensemble.
    if(length(results$pred) == 1) {
      png(file_name, 1600, 1200, res = 350)
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
    x <- (results$pred[[i]]$data$estimate * rate)
    y_test <- (results$pred[[i]]$data$true * rate)
    observed_test <- (results$pred[[i]]$data$observed)
    
    score <- x > cutoff_years
    if(use_hr) {
      sdiff <- tryCatch(survdiff(Surv(y_test, observed_test) ~ score + hr),
                        error = function(e) NULL)
    } else {
      sdiff <- tryCatch(survdiff(Surv(y_test, observed_test) ~ score),
                        error = function(e) NULL)
    }
    
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
    
    pred_labels <- c(paste("<", cutoff_years, "years"), 
                     paste(">=", cutoff_years, "years"))
    pred_survival = pred_labels[1 + (score * 1)]
    df <- data.frame(y = y_test, 
                     observed = observed_test,
                     pred_survival = pred_survival)
    if(use_hr) {
      sfit <- survfit(Surv(y, observed) ~ pred_survival + hr, data = df)
    } else {
      sfit <- survfit(Surv(y, observed) ~ pred_survival, data = df)
    }
    
    # main <- paste("HR and ", title, sep = "")
    
    xname <- result_name
    if(length(grep("clinical", xname, ignore.case = TRUE)) == 1) {
      main <- "null"
    } else if(length(grep("ensemble", xname, ignore.case = TRUE)) == 1) {
      main <- "ensemble"
    } else {
      xname <- gsub("_", " ", xname)
      main <- paste(xname, method_names[i], sep = " - ")
    }
    
    pred_labels <- c(paste0("pred_survival=", pred_labels[1]),
                     paste0("pred_survival=", pred_labels[2]))
    if(use_hr) {
      colors <- matrix(c(paste0(pred_labels[1], ", hr=FALSE"), "blue",
                         paste0(pred_labels[1], ", hr=TRUE "), "orange",
                         paste0(pred_labels[2], ", hr=FALSE"), "black",
                         paste0(pred_labels[2], ", hr=TRUE "), "red"),
                       nrow = 4, byrow = TRUE)
    } else {
      colors <- matrix(c(pred_labels[1], "red",
                         pred_labels[2], "black"),
                       nrow = 2, byrow = TRUE)
    }
    
    if(use_hr) {
      title <- paste0(main, " (", sfit$strata[1], ":", sfit$strata[3], ")")
    } else {
      title <- paste0(main, " (", sfit$strata[1], ":", sfit$strata[2], ")")
    }
    plot(sfit, main = title, 
         #conf.int = TRUE, 
         mark.time = TRUE,
         xmax = 8,
         xlab = "Time (in years)", 
         ylab = y_lab,
         col = colors[colors[, 1] %in% names(sfit$strata), 2])
    if(use_hr) {
      abline(v = 2, col = rgb(0.5, 0.5, 0.5, 0.1), lty = 3)
    }
    #text(x = 1.5, y = 0.05, labels = paste("p =", round(pvalue, 5)))
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
      # Used in poster:
      # plot.km.all(temp, cutoff_years, get_file_name(names(temp)), use_hr = TRUE)
      # plot.km.all(temp, cutoff_years, get_file_name(names(temp)), use_hr = FALSE)
    }
  }
}



