#Please switch "path" with the actual download direction
setwd("path\satellite_data")

# Load datasets
data_files <- c(
  "hstA_converted.csv",
  "hstA_05_converted.csv",
  "graceQ_05.csv",
  "graceA_05.csv",
  "graceQ.csv",
  "hstQ_05.csv",
  "hstQ.csv"
)

rmse_results    <- list()
nlpd_results    <- list()
runtime_results <- list()

m.pred   <- 140
m.hybrid <- 30
m.est    <- 30
g        <- 1 / 10000000
n_splits <- 5  

scale01 <- function(x, min_val, max_val) {
  (x - min_val) / (max_val - min_val)
}


for (file in data_files) {
  cat("\nProcessing dataset:", file, "\n")
  

  df <- read.csv(file, header = TRUE)
  
  df$Time.s. <- df$Time.s. / 10

  X <- as.matrix(df[, c(14, 12, 10)]) 
  y <- as.matrix(df$Cd)               
  t <- as.matrix(df$Time.s.)
  
  kf <- caret::createFolds(y, k = n_splits, list = TRUE, returnTrain = TRUE)
  
  method_rmse    <- list()
  method_nlpd    <- list()
  method_runtime <- list()
  
  for (i in 1:n_splits) {

    train_indices <- unlist(kf[[i]])
    test_indices  <- setdiff(1:nrow(df), train_indices)
    
    trainX <- X[train_indices, ]
    trainY <- y[train_indices]
    testX  <- X[test_indices, ]
    testY  <- y[test_indices]
    trainT <- t[train_indices, ]
    testT  <- t[test_indices, ]
    
    thinning_number <- computeThinningNumber(cbind(trainX, trainY), 30)
    
    ## ---- twinGP ----
    start_time  <- Sys.time()
    fit         <- twingp::twingp(trainX, trainY, testX)
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    preds       <- fit$mu
    vars        <- fit$sigma^2
    rmse        <- sqrt(mean((preds - testY)^2))
    nlpd        <- 0.5 * mean(((testY - preds)^2) / abs(vars) + log(2 * pi * abs(vars)))
    method_rmse$twingp    <- c(method_rmse$twingp, rmse)
    method_nlpd$twingp    <- c(method_nlpd$twingp, nlpd)
    method_runtime$twingp <- c(method_runtime$twingp, elapsed_time)
    
    ## ---- thinned twinGP (full) ----
    start_time  <- Sys.time()
    fit         <- thinned_twingp_full(trainX, trainY, testX, T = thinning_number)
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    preds       <- fit$unweighted$mu
    vars        <- fit$unweighted$sigma^2
    pred_train  <- thinned_twingp_full(trainX, trainY, trainX, T = thinning_number)
    resid       <- trainY - pred_train$unweighted$mu
    pred_g      <- computeLocalFunction(residual = resid,
                                        traindataT = trainT,
                                        testdataT  = testT,
                                        neighbourhood = thinning_number) 
    preds <- preds + pred_g$mean
    vars  <- vars  + pred_g$variance
    rmse  <- sqrt(mean((preds - testY)^2))
    nlpd  <- 0.5 * mean(((testY - preds)^2) / abs(vars) + log(2 * pi * abs(vars)))
    method_rmse$thinned_twingp_full    <- c(method_rmse$thinned_twingp_full, rmse)
    method_nlpd$thinned_twingp_full    <- c(method_nlpd$thinned_twingp_full, nlpd)
    method_runtime$thinned_twingp_full <- c(method_runtime$thinned_twingp_full, elapsed_time)
    
    ## ---- SVecchia ----
    start_time  <- Sys.time()
    fit         <- fit_scaled(trainY, trainX, ms = m.est, max.it = 8,
                              tol.dec = 4, scale = "parms", find.vcf = TRUE)
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    preds_obj   <- predictions_scaled(fit = fit, locs_pred = testX,
                                      nsims = 50, m = m.pred, scale = "parms")
    vars        <- apply(preds_obj$samples, 1, var)
    preds       <- preds_obj$means
    rmse        <- sqrt(mean((preds - testY)^2))
    nlpd        <- 0.5 * mean(((testY - preds)^2) / abs(vars) + log(2 * pi * abs(vars)))
    method_rmse$SVecchia    <- c(method_rmse$SVecchia, rmse)
    method_nlpd$SVecchia    <- c(method_nlpd$SVecchia, nlpd)
    method_runtime$SVecchia <- c(method_runtime$SVecchia, elapsed_time)
    
    ## ---- thinned SV ----
    thinned_bins <- createThinnedBins(trainX, trainY, thinning_number)
    start_time   <- Sys.time()
    fit          <- fit_scaled_thinned(trainY, trainX, scale = "parms",
                                       ms = m.est, max.it = 32,
                                       thinnedBins = thinned_bins,
                                       T = thinning_number, find.vcf = TRUE)
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    preds_obj    <- predictions_scaled_thinned(fit = fit, locs_pred = testX,
                                               nsims = 50, m = m.pred,
                                               scale = "parms")
    vars         <- apply(preds_obj$samples, 1, var)
    preds        <- preds_obj$means
    pred_train   <- predictions_scaled_thinned(fit = fit, locs_pred = trainX,
                                               nsims = 50, m = m.pred,
                                               scale = "parms") 
    resid        <- trainY - pred_train$means
    pred_g       <- computeLocalFunction(residual = resid,
                                         traindataT = trainT,
                                         testdataT  = testT,
                                         neighbourhood = thinning_number) 
    preds <- preds + pred_g$mean
    vars  <- vars  + pred_g$variance
    rmse  <- sqrt(mean((preds - testY)^2))
    nlpd  <- 0.5 * mean(((testY - preds)^2) / abs(vars) + log(2 * pi * abs(vars)))
    method_rmse$thinnedSV    <- c(method_rmse$thinnedSV, rmse)
    method_nlpd$thinnedSV    <- c(method_nlpd$thinnedSV, nlpd)
    method_runtime$thinnedSV <- c(method_runtime$thinnedSV, elapsed_time)
    
    ## ---- SVecchia + t ----
    trainX_t <- cbind(trainT, trainX)
    testX_t  <- cbind(testT,  testX)
    start_time  <- Sys.time()
    fit         <- fit_scaled(trainY, trainX_t, tol.dec = 4, max.it = 8,
                              ms = m.est, scale = "parms", find.vcf = TRUE)
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    preds_obj   <- predictions_scaled(fit = fit, locs_pred = testX_t,
                                      m = m.pred, scale = "parms", nsims = 50)
    vars        <- apply(preds_obj$samples, 1, var)
    preds       <- preds_obj$means
    rmse        <- sqrt(mean((preds - testY)^2))
    nlpd        <- 0.5 * mean(((testY - preds)^2) / abs(vars) + log(2 * pi * abs(vars)))
    method_rmse$SVecchia_t    <- c(method_rmse$SVecchia_t, rmse)
    method_nlpd$SVecchia_t    <- c(method_nlpd$SVecchia_t, nlpd)
    method_runtime$SVecchia_t <- c(method_runtime$SVecchia_t, elapsed_time)
    
    ## ---- global-local laGP ----
    start_time  <- Sys.time()
    preds_obj   <- global_local_laGP1(trainX, trainY, testX, m.hybrid)
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    preds       <- preds_obj$mean
    vars        <- preds_obj$var
    rmse        <- sqrt(mean((preds - testY)^2))
    nlpd        <- 0.5 * mean(((testY - preds)^2) / abs(vars) + log(2 * pi * abs(vars)))
    method_rmse$global_local_laGP    <- c(method_rmse$global_local_laGP, rmse)
    method_nlpd$global_local_laGP    <- c(method_nlpd$global_local_laGP, nlpd)
    method_runtime$global_local_laGP <- c(method_runtime$global_local_laGP, elapsed_time)
    
    ## ---- thinned laGP ----
    start_time  <- Sys.time()
    preds_obj   <- thinned_lagp(trainX, trainY, testX)
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    preds       <- preds_obj$mu
    vars        <- (preds_obj$sigma)^2
    rmse        <- sqrt(mean((preds - testY)^2))
    nlpd        <- 0.5 * mean(((testY - preds)^2) / abs(vars) + log(2 * pi * abs(vars)))
    method_rmse$thinned_lagp    <- c(method_rmse$thinned_lagp, rmse)
    method_nlpd$thinned_lagp    <- c(method_nlpd$thinned_lagp, nlpd)
    method_runtime$thinned_lagp <- c(method_runtime$thinned_lagp, elapsed_time)
    
    

    x_min <- apply(trainX, 2, min)
    x_max <- apply(trainX, 2, max)
    
    trainX_scaled <- sweep(trainX, 2, x_min, FUN = "-")
    trainX_scaled <- sweep(trainX_scaled, 2, x_max - x_min, FUN = "/")
    
    testX_scaled  <- sweep(testX, 2, x_min, FUN = "-")
    testX_scaled  <- sweep(testX_scaled, 2, x_max - x_min, FUN = "/")
    
    ## ---- tempGP ----
    start_time <- Sys.time()
    fit <- tempGP(as.matrix(trainX_scaled), as.vector(trainY), trainT,
                  fast_computation   = FALSE,
                  max_thinning_number = 30L,
                  limit_memory        = 5000L)
    
    preds <- predict(fit, as.matrix(testX_scaled), testT)
    
    elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
    rmse <- sqrt(mean((preds - testY)^2))
    nlpd <- 0   
    method_rmse$tempGP    <- c(method_rmse$tempGP, rmse)
    method_nlpd$tempGP    <- c(method_nlpd$tempGP, nlpd)
    method_runtime$tempGP <- c(method_runtime$tempGP, elapsed_time)
    
  } 
  
  rmse_results[[file]]    <- method_rmse
  nlpd_results[[file]]    <- method_nlpd
  runtime_results[[file]] <- method_runtime
}

# Save intermediate results
#write.xlsx(rmse_results,    file = "RMSE_Results_12_21.xlsx")
#write.xlsx(nlpd_results,    file = "NLPD_Results_12_21.xlsx")
#write.xlsx(runtime_results, file = "Runtime_Results_12_21.xlsx")

#print("Results saved to Excel files.")


methods_order <- c(
  "global_local_laGP",   # laGP
  "twingp",              # twinGP
  "SVecchia",            # SV(x)
  "SVecchia_t",          # SV(x,t)
  "thinned_lagp",        # thinned laGP
  "thinned_twingp_full", # thinned twinGP
  "thinnedSV",           # thinned SV
  "tempGP"
)

method_labels <- c(
  global_local_laGP   = "laGP",
  twingp              = "twinGP",
  SVecchia            = "SV(x)",
  SVecchia_t          = "SV(x,t)",
  thinned_lagp        = "thinned laGP",
  thinned_twingp_full = "thinned twinGP",
  thinnedSV           = "thinned SV",
  tempGP              = "tempGP"
)

datasets_order <- c(
  "hstA_converted.csv",
  "hstA_05_converted.csv",
  "hstQ.csv",
  "hstQ_05.csv",
  "graceA_05.csv",
  "graceQ.csv",
  "graceQ_05.csv"
)

dataset_labels <- c(
  "hstA_converted.csv"    = "hstA",
  "hstA_05_converted.csv" = "hstA_05",
  "hstQ.csv"              = "hstQ",
  "hstQ_05.csv"           = "hstQ_05",
  "graceA_05.csv"         = "graceA_05",
  "graceQ.csv"            = "graceQ",
  "graceQ_05.csv"         = "graceQ_05"
)


rmse_mat <- sapply(datasets_order, function(dset) {
  sapply(methods_order, function(meth) {
    mean(rmse_results[[dset]][[meth]])
  })
})

rmse_mat <- t(rmse_mat)
rmse_mat <- rmse_mat[, methods_order, drop = FALSE]
rmse_mat <- t(rmse_mat)

rownames(rmse_mat) <- method_labels[methods_order]
colnames(rmse_mat) <- dataset_labels[datasets_order]

rmse_mat <- rmse_mat * 1000
rmse_table <- cbind(rmse_mat, Average = rowMeans(rmse_mat))
table3_disp <- round(rmse_table, 2)


nlpd_mat <- sapply(datasets_order, function(dset) {
  sapply(methods_order, function(meth) {
    mean(nlpd_results[[dset]][[meth]])
  })
})

nlpd_mat <- t(nlpd_mat)
nlpd_mat <- nlpd_mat[, methods_order, drop = FALSE]
nlpd_mat <- t(nlpd_mat)

rownames(nlpd_mat) <- method_labels[methods_order]
colnames(nlpd_mat) <- dataset_labels[datasets_order]

nlpd_table <- cbind(nlpd_mat, Average = rowMeans(nlpd_mat))
table4_disp <- round(nlpd_table, 2)

runtime_sec <- sapply(methods_order, function(meth) {
  vals <- unlist(lapply(runtime_results, function(d) {
    as.numeric(d[[meth]], units = "secs")
  }))
  mean(vals)
})

runtime_table <- t(as.matrix(runtime_sec))
colnames(runtime_table) <- method_labels[methods_order]
rownames(runtime_table) <- "Runtime (sec)"

table5_disp <- round(runtime_table, 0)

table3_disp
table4_disp
table5_disp
