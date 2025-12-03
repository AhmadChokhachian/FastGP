#Please switch "path" with the actual download direction
setwd("path/Dataset5")

datasets <- c(
  "Inland Wind Farm Dataset1(WT1)", 
  "Inland Wind Farm Dataset1(WT2)",
  "Inland Wind Farm Dataset1(WT3)", 
  "Inland Wind Farm Dataset1(WT4)",
  "Offshore Wind Farm Dataset1(WT5)",
  "Offshore Wind Farm Dataset1(WT6)"
)

# Settings
m.pred   <- 140
m.hybrid <- 30
m.est    <- 30
g        <- 1 / 10000000
n_splits <- 5  


rmse_results    <- list()
nlpd_results    <- list()
runtime_results <- list()

for (dataset in datasets) {
  df <- read.csv(paste0(dataset, ".csv"), header = TRUE, sep = ",")
  set.seed(1)
  
  # Add wind direction sin and cos columns
  df$wind_direction_sin <- sin(df$D)
  df$wind_direction_cos <- cos(df$D)
  
  if (dataset %in% c("Offshore Wind Farm Dataset1(WT5)",
                     "Offshore Wind Farm Dataset1(WT6)")) {
    Ycol <- 9
    Tcol <- 1
    xCol <- c(2, 3, 4, 5, 6, 7, 8, 10, 11)
  } else {
    Ycol <- 7
    Tcol <- 1
    xCol <- c(2, 3, 4, 5, 6, 8, 9)
  }
  
  
  kf <- caret::createFolds(df[, Ycol], k = n_splits,
                           list = TRUE, returnTrain = TRUE)
  
  method_rmse    <- list()
  method_nlpd    <- list()
  method_runtime <- list()
  
  for (i in 1:n_splits) {
    
    train_indices <- unlist(kf[[i]])
    test_indices  <- setdiff(1:nrow(df), train_indices)
    
    
    trainX <- as.matrix(df[train_indices, xCol])
    trainY <- as.numeric(df[train_indices, Ycol])
    testX  <- as.matrix(df[test_indices,  xCol])
    testY  <- as.numeric(df[test_indices, Ycol])
    trainT <- as.matrix(df[train_indices, Tcol])
    testT  <- as.matrix(df[test_indices,  Tcol])
    
    ### Methods ###
    
    ## twinGP
    start_time_fit <- system.time({
      fit <- twingp::twingp(trainX, trainY, testX)
    })["elapsed"]
    start_time_fit <- as.numeric(start_time_fit)
    preds <- fit$mu
    vars  <- fit$sigma^2
    rmse  <- sqrt(mean((preds - testY)^2))
    nlpd  <- 0.5 * mean(((testY - preds)^2) / abs(vars) +
                          log(2 * pi * abs(vars)))
    method_runtime$twingp <- c(method_runtime$twingp, start_time_fit)
    method_rmse$twingp    <- c(method_rmse$twingp, rmse)
    method_nlpd$twingp    <- c(method_nlpd$twingp, nlpd)
    
    ## global_local_laGP
    start_time_fit <- system.time({
      preds_obj <- global_local_laGP1(trainX, trainY, testX, m.hybrid)
    })["elapsed"]
    start_time_fit <- as.numeric(start_time_fit)
    preds <- preds_obj$mean
    vars  <- preds_obj$var
    rmse  <- sqrt(mean((preds - testY)^2))
    nlpd  <- 0.5 * mean(((testY - preds)^2) / abs(vars) +
                          log(2 * pi * abs(vars)))
    method_runtime$global_local_laGP <- c(method_runtime$global_local_laGP,
                                          start_time_fit)
    method_rmse$global_local_laGP <- c(method_rmse$global_local_laGP, rmse)
    method_nlpd$global_local_laGP <- c(method_nlpd$global_local_laGP, nlpd)
    
    ## SVecchia
    start_time_fit <- system.time({
      fit <- fit_scaled(trainY, trainX, ms = m.est,
                        max.it = 8, tol.dec = 4,
                        scale = "parms", find.vcf = TRUE)
      preds <- predictions_scaled(fit = fit, locs_pred = testX,
                                  nsims = 50, m = m.pred,
                                  scale = "parms")
    })["elapsed"]
    start_time_fit <- as.numeric(start_time_fit)
    vars  <- apply(preds$samples, 1, var)
    preds <- preds$means
    rmse  <- sqrt(mean((preds - testY)^2))
    nlpd  <- 0.5 * mean(((testY - preds)^2) / abs(vars) +
                          log(2 * pi * abs(vars)))
    method_runtime$SVecchia <- c(method_runtime$SVecchia, start_time_fit)
    method_rmse$SVecchia    <- c(method_rmse$SVecchia, rmse)
    method_nlpd$SVecchia    <- c(method_nlpd$SVecchia, nlpd)
    
    ## thinnedSV
    thinning_number <- computeThinningNumber(trainX, 50)
    thinned_bins    <- createThinnedBins(trainX, trainY, thinning_number)
    start_time_fit  <- system.time({
      fit <- fit_scaled_thinned(trainY, trainX, scale = "parms",
                                ms = m.est,
                                thinnedBins = thinned_bins,
                                T = thinning_number,
                                find.vcf = TRUE)
      preds <- predictions_scaled_thinned(fit = fit, locs_pred = testX,
                                          nsims = 50, m = m.pred,
                                          scale = "parms")
    })["elapsed"]
    start_time_fit <- as.numeric(start_time_fit)
    vars  <- apply(preds$samples, 1, var)
    preds <- preds$means
    pred_train <- predictions_scaled_thinned(fit = fit, locs_pred = trainX,
                                             nsims = 50, m = m.pred,
                                             scale = "parms")
    resid <- trainY - pred_train$means
    pred_g <- computeLocalFunction(residual   = resid,
                                   traindataT = trainT,
                                   testdataT  = testT,
                                   neighbourhood = thinning_number)
    preds <- preds + pred_g$mean
    vars  <- vars  + pred_g$variance
    rmse  <- sqrt(mean((preds - testY)^2))
    nlpd  <- 0.5 * mean(((testY - preds)^2) / abs(vars) +
                          log(2 * pi * abs(vars)))
    method_runtime$thinnedSV <- c(method_runtime$thinnedSV, start_time_fit)
    method_rmse$thinnedSV    <- c(method_rmse$thinnedSV, rmse)
    method_nlpd$thinnedSV    <- c(method_nlpd$thinnedSV, nlpd)
    
    ## SVecchia + t
    trainX_t <- cbind(trainT, trainX)
    testX_t  <- cbind(testT,  testX)
    start_time_fit <- system.time({
      fit <- fit_scaled(trainY, trainX_t, tol.dec = 4, max.it = 8,
                        ms = m.est, scale = "parms", find.vcf = TRUE)
      preds <- predictions_scaled(fit = fit, locs_pred = testX_t,
                                  m = m.pred, scale = "parms",
                                  nsims = 50)
    })["elapsed"]
    start_time_fit <- as.numeric(start_time_fit)
    vars  <- apply(preds$samples, 1, var)
    preds <- preds$means
    rmse  <- sqrt(mean((preds - testY)^2))
    nlpd  <- 0.5 * mean(((testY - preds)^2) / abs(vars) +
                          log(2 * pi * abs(vars)))
    method_runtime$SVecchia_t <- c(method_runtime$SVecchia_t, start_time_fit)
    method_rmse$SVecchia_t    <- c(method_rmse$SVecchia_t, rmse)
    method_nlpd$SVecchia_t    <- c(method_nlpd$SVecchia_t, nlpd)
    
    ## tempGP (no NLPD in final table)
    start_time_fit <- system.time({
      fit <- DSWE::tempGP(as.matrix(trainX), as.vector(trainY), trainT,
                          fast_computation    = FALSE,
                          max_thinning_number = 30L)
      preds <- predict(fit, as.matrix(testX), testT)
    })["elapsed"]
    start_time_fit <- as.numeric(start_time_fit)
    rmse <- sqrt(mean((preds - testY)^2))
    # NLPD is not meaningful here -> record NA
    nlpd <- NA_real_
    method_runtime$tempgp <- c(method_runtime$tempgp, start_time_fit)
    method_rmse$tempgp    <- c(method_rmse$tempgp, rmse)
    method_nlpd$tempgp    <- c(method_nlpd$tempgp, nlpd)
    
    ## thinned twinGP
    thinning_number <- computeThinningNumber(trainX, 20)
    start_time_fit  <- system.time({
      fit <- thinned_twingp_full(trainX, trainY, testX, T = thinning_number)
    })["elapsed"]
    start_time_fit <- as.numeric(start_time_fit)
    pred_train    <- thinned_twingp_full(trainX, trainY, trainX,
                                         T = thinning_number)
    pred_train_mu <- pred_train$unweighted$mu
    resid         <- trainY - pred_train_mu
    pred_g <- computeLocalFunction(residual   = resid,
                                   traindataT = trainT,
                                   testdataT  = testT,
                                   neighbourhood = thinning_number)
    preds <- fit$unweighted$mu + pred_g$mean
    vars  <- fit$unweighted$sigma^2 + pred_g$variance
    rmse  <- sqrt(mean((preds - testY)^2))
    nlpd  <- 0.5 * mean(((testY - preds)^2) / abs(vars) +
                          log(2 * pi * abs(vars)))
    method_runtime$thinned_twingp_full <- c(method_runtime$thinned_twingp_full,
                                            start_time_fit)
    method_rmse$thinned_twingp_full <- c(method_rmse$thinned_twingp_full, rmse)
    method_nlpd$thinned_twingp_full <- c(method_nlpd$thinned_twingp_full, nlpd)
    
    ## thinned laGP
    thinning_number <- computeThinningNumber(trainX, 20)
    start_time_fit  <- system.time({
      fit <- thinned_lagp(trainX, trainY, testX, T = thinning_number)
    })["elapsed"]
    start_time_fit <- as.numeric(start_time_fit)
    pred_train    <- thinned_lagp(trainX, trainY, trainX,
                                  T = thinning_number)
    pred_train_mu <- pred_train$mu
    resid         <- trainY - pred_train_mu
    pred_g <- computeLocalFunction(residual   = resid,
                                   traindataT = trainT,
                                   testdataT  = testT,
                                   neighbourhood = thinning_number)
    preds <- fit$mu + pred_g$mean
    vars  <- fit$sigma^2 + pred_g$variance
    rmse  <- sqrt(mean((preds - testY)^2))
    nlpd  <- 0.5 * mean(((testY - preds)^2) / abs(vars) +
                          log(2 * pi * abs(vars)))
    method_runtime$thinned_lagp_full <- c(method_runtime$thinned_lagp_full,
                                          start_time_fit)
    method_rmse$thinned_lagp_full <- c(method_rmse$thinned_lagp_full, rmse)
    method_nlpd$thinned_lagp_full <- c(method_nlpd$thinned_lagp_full, nlpd)
    
  } # end fold loop
  
  # Store results for this dataset
  rmse_results[[dataset]]    <- method_rmse
  nlpd_results[[dataset]]    <- method_nlpd
  runtime_results[[dataset]] <- method_runtime
}

###########################################
## Build Tables 6, 7, and 9 for Dataset 5
###########################################

# Mapping dataset names to WT labels
wt_map <- c(
  "Inland Wind Farm Dataset1(WT1)"   = "WT1",
  "Inland Wind Farm Dataset1(WT2)"   = "WT2",
  "Inland Wind Farm Dataset1(WT3)"   = "WT3",
  "Inland Wind Farm Dataset1(WT4)"   = "WT4",
  "Offshore Wind Farm Dataset1(WT5)" = "WT5",
  "Offshore Wind Farm Dataset1(WT6)" = "WT6"
)
wt_order <- c("WT1","WT2","WT3","WT4","WT5","WT6")

# Internal method names and pretty labels
methods_internal <- c(
  "global_local_laGP",
  "twingp",
  "SVecchia",
  "SVecchia_t",
  "thinned_lagp_full",
  "thinned_twingp_full",
  "thinnedSV",
  "tempgp"
)

method_labels <- c(
  global_local_laGP   = "laGP",
  twingp              = "twinGP",
  SVecchia            = "SV(x)",
  SVecchia_t          = "SV(x,t)",
  thinned_lagp_full   = "thinned laGP",
  thinned_twingp_full = "thinned twinGP",
  thinnedSV           = "thinned SV",
  tempgp              = "tempGP"
)


rmse_mat <- sapply(datasets, function(dset) {
  sapply(methods_internal, function(m) {
    mean(rmse_results[[dset]][[m]])
  })
})

colnames(rmse_mat) <- wt_map[colnames(rmse_mat)]
rmse_mat <- rmse_mat[, wt_order, drop = FALSE]
rmse_mat <- t(rmse_mat)  
rmse_mat <- t(rmse_mat)  

# Average RMSE across WT1–WT6
rmse_avg <- rowMeans(rmse_mat)

table6_df <- data.frame(
  Method = method_labels[methods_internal],
  RMSE_WT1 = rmse_mat[,"WT1"],
  RMSE_WT2 = rmse_mat[,"WT2"],
  RMSE_WT3 = rmse_mat[,"WT3"],
  RMSE_WT4 = rmse_mat[,"WT4"],
  RMSE_WT5 = rmse_mat[,"WT5"],
  RMSE_WT6 = rmse_mat[,"WT6"],
  Average  = rmse_avg,
  row.names = NULL
)


nlpd_mat <- sapply(datasets, function(dset) {
  sapply(methods_internal, function(m) {
    mean(nlpd_results[[dset]][[m]])
  })
})
colnames(nlpd_mat) <- wt_map[colnames(nlpd_mat)]
nlpd_mat <- nlpd_mat[, wt_order, drop = FALSE]
nlpd_mat <- t(nlpd_mat)
nlpd_mat <- t(nlpd_mat)


nlpd_mat["tempgp", ] <- NA_real_

nlpd_avg <- rowMeans(nlpd_mat, na.rm = FALSE) 

table7_df <- data.frame(
  Method = method_labels[methods_internal],
  NLPD_WT1 = nlpd_mat[,"WT1"],
  NLPD_WT2 = nlpd_mat[,"WT2"],
  NLPD_WT3 = nlpd_mat[,"WT3"],
  NLPD_WT4 = nlpd_mat[,"WT4"],
  NLPD_WT5 = nlpd_mat[,"WT5"],
  NLPD_WT6 = nlpd_mat[,"WT6"],
  Average  = nlpd_avg,
  row.names = NULL
)

# ----- Runtime – Table 9 (Dataset 5 row, minutes) -----

runtime_sec <- sapply(methods_internal, function(m) {
  
  all_vals <- unlist(lapply(runtime_results, function(one_dataset) {
    one_dataset[[m]]
  }))
  mean(all_vals)
})

runtime_min <- runtime_sec / 60

table9_df <- data.frame(
  Dataset = "Dataset 5",
  t(as.data.frame(as.list(runtime_min))),
  row.names = NULL
)


colnames(table9_df) <- c("Dataset", method_labels[methods_internal])


table6_disp <- table6_df
table6_disp[,-1] <- round(table6_disp[,-1], 2)

table7_disp <- table7_df
table7_disp[,-1] <- round(table7_disp[,-1], 2)

table9_disp <- table9_df
table9_disp[,-1] <- round(table9_disp[,-1], 1)

cat("\n=== Table 6 (RMSE, DSWE Dataset 5) ===\n")
print(table6_disp)

cat("\n\n=== Table 7 (NLPD, DSWE Dataset 5) ===\n")
print(table7_disp)

cat("\n\n=== Table 9 (Runtime, Dataset 5, minutes) ===\n")
print(table9_disp)


#write.xlsx(table6_disp, "Table6_DSWE5_RMSE.xlsx", row.names = FALSE)
#write.xlsx(table7_disp, "Table7_DSWE5_NLPD.xlsx", row.names = FALSE)
#write.xlsx(table9_disp, "Table9_DSWE5_Runtime.xlsx", row.names = FALSE)
