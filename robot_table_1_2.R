setwd("robot_data")
 
run_experiment <- function(version, methods_to_run, reps_to_run,
                           m.est = 30, m.pred = 140, n.est = 3000,
                           data_dir = "data", save_dir = "results") {
  
  library(tictoc)
  dir.create(save_dir, showWarnings = FALSE)
  
  all_methods <- c("SVecchia_t", "SVecchia", "Hybrid_laGP", "naive_twingp",
                   "SVecchia_thinned", "thinned_twingp_full",
                   "thinned_lagp", "tempGP")
  
  method_indices <- match(methods_to_run, all_methods)
  
  rmse <- matrix(NA, nrow = length(all_methods), ncol = max(reps_to_run))
  nlpd <- matrix(NA, nrow = length(all_methods), ncol = max(reps_to_run))
  time <- matrix(NA, nrow = length(all_methods), ncol = max(reps_to_run))
  rownames(rmse) <- all_methods
  rownames(nlpd) <- all_methods
  rownames(time) <- all_methods
  
  source('https://raw.githubusercontent.com/katzfuss-group/scaledVecchia/master/test_functions.R')
  fun <- testfuns[[2]]$fun
  
  for (rep in reps_to_run) {
    cat(sprintf("Running version=%s, rep=%d...\n", version, rep))
    
    # Load data
    train_path <- file.path(data_dir, version, sprintf("rep_%d_train.csv", rep))
    test_path <- file.path(data_dir, version, sprintf("rep_%d_test.csv", rep))
    
    train_df <- read.csv(train_path)
    test_df <- read.csv(test_path)
    test_df=test_df[1:10000,]
    inputs.train <- as.matrix(train_df[, 1:8])
    print(sum(is.na(inputs.train)))
    y.train <- train_df$y
    inputs.test <- as.matrix(test_df[, 1:8])
    y.test <- test_df$y
    
    # Compute T if needed
    if (any(methods_to_run %in% c("SVecchia_thinned", "thinned_twingp", "thinned_twingp_full", "thinned_lagp", "thinned_lagp_full"))) {
      T <- computeThinningNumber(cbind(inputs.train,y.train), 30)
      print(T)
    }
    
    for (method in methods_to_run) {
      idx <- match(method, all_methods)
      t1 <- tic()
      
      # ==== Implementations ====
      if (method == "SVecchia_t") {
        idx.train <- seq_len(nrow(inputs.train))
        idx.test <- seq_len(nrow(inputs.test))
        inputs.train_t <- cbind(idx.train, inputs.train)
        inputs.test_t <- cbind(idx.test, inputs.test)
        
        fit <- fit_scaled(y.train, inputs.train_t, ms = m.est, scale = "parms", trend = 'zero',
                          #n.est = n.est,
                          print.level = 0, find.vcf = TRUE, nu = 1.5, nug = 0,max.it = 8)
        preds <- predictions_scaled(fit, locs_pred = inputs.test_t, m = m.pred, scale = "parms", nsims = 50)
        yhat <- preds$means
        vars <- apply(preds$samples, 1, var)
        
      } else if (method == "SVecchia") {
        fit <- fit_scaled(y.train, inputs.train, ms = m.est, scale = "parms", trend = 'zero',
                          #n.est = n.est, 
                          print.level = 0, find.vcf = TRUE, nu = 1.5, nug = 0,max.it = 8)
        preds <- predictions_scaled(fit, locs_pred = inputs.test, m = m.pred, scale = "parms", nsims = 50)
        yhat <- preds$means
        vars <- apply(preds$samples, 1, var)
        
      } else if (method == "Hybrid_laGP") {
        scale_01 <- function(x, min_val, max_val) {
          (x - min_val) / (max_val - min_val)
        }
        
        # Compute min and max from training inputs
        min_vals <- apply(inputs.train, 2, min)
        max_vals <- apply(inputs.train, 2, max)
        
        # Scale training and test inputs
        inputs.train.scaled <- sweep(sweep(inputs.train, 2, min_vals, "-"), 2, max_vals - min_vals, "/")
        inputs.test.scaled  <- sweep(sweep(inputs.test,  2, min_vals, "-"), 2, max_vals - min_vals, "/")
        
        pred <- global_local_laGP1(inputs.train.scaled, y.train, inputs.test.scaled, m.est)
        yhat <- pred$mean
        vars <- pred$var
        
      } else if (method == "naive_twingp") {
        pred <- twingp(inputs.train, y.train, inputs.test, g_num = m.pred)
        yhat <- pred$mu
        vars <- pred$sigma^2
        
      } else if (method == "SVecchia_thinned") {
        bins <- createThinnedBins(inputs.train, y.train, thinningNumber = T)
        fit <- fit_scaled_thinned(y.train, inputs.train, ms = m.est, trend = 'zero', scale = "parms", T = T,
                                  thinnedBins = bins,max.it = 32)
        pred <- predictions_scaled_thinned(fit, locs_pred = inputs.test, m = m.pred, scale = "parms", nsims = 50)
        yhat <- pred$means
        vars <- apply(pred$samples, 1, var)
        
      } else if (method == "thinned_twingp") {
        pred <- thinned_twingp(inputs.train, y.train, inputs.test, T = T)
        yhat <- pred$mu
        vars <- pred$sigma^2
        
      } else if (method == "thinned_twingp_full") {
        pred <- thinned_twingp_full(inputs.train, y.train, inputs.test, T = T)
        yhat <- pred$unweighted$mu
        vars <- pred$unweighted$sigma^2
        
      } else if (method == "thinned_lagp") {
        pred <- thinned_lagp(inputs.train, y.train, inputs.test, T = T)
        yhat <- pred$mu
        vars <- pred$sigma^2
        
      } else if (method == "thinned_lagp_full") {
        pred <- thinned_lagp_full(inputs.train, y.train, inputs.test, T = T)
        yhat <- pred$unweighted$mu
        vars <- pred$unweighted$sigma^2
        
      } else if (method == "tempGP") {
        # Manually scale inputs to [0, 1] range
        scale_01 <- function(x, min_val, max_val) {
          (x - min_val) / (max_val - min_val)
        }
        
        # Compute min and max from training inputs
        min_vals <- apply(inputs.train, 2, min)
        max_vals <- apply(inputs.train, 2, max)
        
        # Scale training and test inputs
        inputs.train.scaled <- sweep(sweep(inputs.train, 2, min_vals, "-"), 2, max_vals - min_vals, "/")
        inputs.test.scaled  <- sweep(sweep(inputs.test,  2, min_vals, "-"), 2, max_vals - min_vals, "/")
        
        # Fit the tempGP model on scaled inputs
        fit <- DSWE::tempGP(inputs.train.scaled, y.train,
                            max_thinning_number = 30L,
                            fast_computation = FALSE, limit_memory = 5000L)
        
        # Predict on scaled test inputs
        yhat <- predict(fit, inputs.test.scaled)
        
        vars <- rep(NA, length(yhat))  # no uncertainty
      }
      
      temp <- toc(quiet = TRUE)
      time[idx, rep] <- temp$toc - temp$tic
      rmse[idx, rep] <- sqrt(mean((yhat - y.test)^2))
      print(rmse)
      print(time)
      print(rmse)
      if (!all(is.na(vars))) {
        nlpd[idx, rep] <- 0.5 * mean(((y.test - yhat)^2) / vars + log(2 * pi * vars))
      } else {
        nlpd[idx, rep] <- NA  # Not available (e.g., for tempGP)
      }
    }
  }
  
  # Save outputs
  rep_label <- paste0("reps_", min(reps_to_run), "-", max(reps_to_run))
  method_label <- if (length(methods_to_run) == length(all_methods)) "allMethods" else paste(methods_to_run, collapse = "_")
  
  
  return(list(
    rmse = rmse,
    nlpd = nlpd,
    time = time,
    methods = all_methods
  ))
  
#  write.csv(rmse, file = file.path(save_dir, sprintf("rmse_%s_%s_%s.csv", version, method_label, rep_label)))
#  write.csv(nlpd, file = file.path(save_dir, sprintf("nlpd_%s_%s_%s.csv", version, method_label, rep_label)))
#  write.csv(time, file = file.path(save_dir, sprintf("time_%s_%s_%s.csv", version, method_label, rep_label)))
}

###low
out_low  <- run_experiment(
  version = "version5",
  methods_to_run = c("SVecchia_t", "SVecchia", "Hybrid_laGP", "naive_twingp",
                                     "SVecchia_thinned", "thinned_twingp_full",
                                     "thinned_lagp", "tempGP"),
  reps_to_run = 1:10
)
###medium
out_med  <- run_experiment(
  version = "version2",
  methods_to_run = c("SVecchia_t", "SVecchia", "Hybrid_laGP", "naive_twingp",
                    "SVecchia_thinned", "thinned_twingp_full",
                    "thinned_lagp", "tempGP"),
  reps_to_run = 1:10
)
###high
out_high  <- run_experiment(
  version = "version7",
  methods_to_run = c("SVecchia_t", "SVecchia", "Hybrid_laGP", "naive_twingp",
                    "SVecchia_thinned", "thinned_twingp_full",
                    "thinned_lagp", "tempGP"),
  reps_to_run = 1:10
)


row_mean <- function(x) apply(x, 1, mean, na.rm = TRUE)
row_sd   <- function(x) apply(x, 1, sd,   na.rm = TRUE)

rmse_mean_low  <- row_mean(out_low$rmse)
rmse_sd_low    <- row_sd(out_low$rmse)

rmse_mean_med  <- row_mean(out_med$rmse)
rmse_sd_med    <- row_sd(out_med$rmse)

rmse_mean_high <- row_mean(out_high$rmse)
rmse_sd_high   <- row_sd(out_high$rmse)

nlpd_mean_low  <- row_mean(out_low$nlpd)
nlpd_sd_low    <- row_sd(out_low$nlpd)

nlpd_mean_med  <- row_mean(out_med$nlpd)
nlpd_sd_med    <- row_sd(out_med$nlpd)

nlpd_mean_high <- row_mean(out_high$nlpd)
nlpd_sd_high   <- row_sd(out_high$nlpd)


fmt <- function(m,s) sprintf("%.3f ± %.3f", m, s)

rmse_low_col  <- fmt(rmse_mean_low,  rmse_sd_low)
rmse_med_col  <- fmt(rmse_mean_med,  rmse_sd_med)
rmse_high_col <- fmt(rmse_mean_high, rmse_sd_high)

nlpd_low_col  <- fmt(nlpd_mean_low,  nlpd_sd_low)
nlpd_med_col  <- fmt(nlpd_mean_med,  nlpd_sd_med)
nlpd_high_col <- fmt(nlpd_mean_high, nlpd_sd_high)

time_mean_low  <- row_mean(out_low$time)
time_mean_med  <- row_mean(out_med$time)
time_mean_high <- row_mean(out_high$time)

runtime_avg <- (time_mean_low + time_mean_med + time_mean_high) / 3


method_labels <- c(
  "Hybrid_laGP"         = "laGP",
  "naive_twingp"        = "twinGP",
  "SVecchia"            = "SV(x)",
  "SVecchia_t"          = "SV(x,t)",
  "thinned_lagp"        = "thinned laGP",
  "thinned_twingp_full" = "thinned twinGP",
  "SVecchia_thinned"    = "thinned SV",
  "tempGP"              = "tempGP"
)

order_methods <- c("Hybrid_laGP","naive_twingp","SVecchia","SVecchia_t",
                   "thinned_lagp","thinned_twingp_full",
                   "SVecchia_thinned","tempGP")

Table1 <- data.frame(
  Method     = method_labels[order_methods],
  RMSE_low   = rmse_low_col[order_methods],
  RMSE_med   = rmse_med_col[order_methods],
  RMSE_high  = rmse_high_col[order_methods],
  NLPD_low   = nlpd_low_col[order_methods],
  NLPD_med   = nlpd_med_col[order_methods],
  NLPD_high  = nlpd_high_col[order_methods]
)
table1_disp
Table2 <- data.frame(
  Method = method_labels[order_methods],
  Runtime_sec = round(runtime_avg[order_methods], 0)
)
table2_disp
