setwd("/Dataset6")

datasets <- c("Inland Wind Farm Dataset2(WT1)",
              "Inland Wind Farm Dataset2(WT2)", 
              "Offshore Wind Farm Dataset2(WT3)", 
              "Offshore Wind Farm Dataset2(WT4)"
)

results_df <- data.frame(
  Dataset     = character(),
  Method      = character(),
  RMSE_T2     = numeric(),
  Runtime_Fit = numeric(),
  Runtime_T2  = numeric(),
  NLPD_T2     = numeric(),
  stringsAsFactors = FALSE
)

m.pred   <- 140
m.hybrid <- 30
m.est    <- 30

# choose whichever methods you want to run
methods <- c(
  "H-laGP",
  "twingp",
  "SVecchia",
  "SVecchia+t",
  "thinned_lagp",
  "thinned_twingp_full",
  "thinnedSV",
  "tempgp"
)

for (dataset in datasets) {
  df <- read.csv(paste0(dataset, ".csv"), header = TRUE, sep = ",")
  set.seed(1)
  
  df$wind_direction_sin <- sin(df$D)
  df$wind_direction_cos <- cos(df$D)
  
  xCol <- c(2, 3, 4, 5, 6, 7, 9, 10) 
  yCol <- 8                          
  
  # Split data into T1, T2
  if (dataset %in% c("Offshore Wind Farm Dataset2(WT3)",
                     "Offshore Wind Farm Dataset2(WT4)")) {
    t1 <- subset(df, year(time) %in% c(2007, 2008))
    t2 <- subset(df, year(time) == 2009)
  } else {
    t1 <- subset(df, year(time) %in% c(2008, 2009))
    t2 <- subset(df, year(time) == 2010)
  }
  
  trainX    <- as.matrix(t1[, xCol])
  trainY    <- as.numeric(t1[, yCol])
  test_T2_X <- as.matrix(t2[, xCol])
  
  for (method in methods) {
    cat("Dataset:", dataset, "Method:", method, "\n")
    

    fit_runtime <- system.time({
      if (method == "SVecchia") {
        scale <- "parms"
        fit <- fit_scaled(trainY, trainX,
                          ms = m.est, scale = scale,
                          trend = "linear",
                          print.level = 2,
                          find.vcf = TRUE)
      } else if (method == "H-laGP") {
        fit <- list() 
      } else if (method == "thinnedSV") {
        Tthin <- computeThinningNumber(trainX, 20)
        fit <- fit_scaled_thinned(trainY, trainX,
                                  ms = m.est,
                                  thinnedBins = createThinnedBins(trainX, trainY, Tthin),
                                  T = Tthin, find.vcf = TRUE)
      } else if (method == "twingp") {
        fit <- list(trainX = trainX, trainY = trainY)
      } else if (method == "thinned_twingp_full") {
        fit <- list(trainX = trainX, trainY = trainY)
      } else if (method == "thinned_lagp") {
        fit <- list(trainX = trainX, trainY = trainY)
      } else if (method == "tempgp") {
        fit <- list(trainX = trainX, trainY = trainY)
      } else if (method == "SVecchia+t") {
        scale <- "parms"
        idx_train   <- 1:nrow(trainX)
        idx_test_T2 <- (nrow(trainX) + 1):(nrow(trainX) + nrow(test_T2_X))
        trainX_t    <- cbind(idx_train,   trainX)
        test_T2_X_t <- cbind(idx_test_T2, test_T2_X)
        fit <- fit_scaled(trainY, trainX_t,
                          ms = m.est, scale = scale,
                          trend = "linear",
                          find.vcf = TRUE)
      }
    })["elapsed"]
    fit_runtime <- as.numeric(fit_runtime) 
    
    runtime_T2 <- system.time({
      if (method == "SVecchia") {
        preds_T2 <- predictions_scaled(fit = fit,
                                       locs_pred = test_T2_X,
                                       m = m.pred,
                                       scale = scale,
                                       nsims = 50)
        means_T2 <- preds_T2$means
        vars_T2  <- apply(preds_T2$samples, 1, var)
      } else if (method == "H-laGP") {
        preds_T2 <- global_local_laGP1(trainX, trainY, test_T2_X, m.hybrid)
        means_T2 <- preds_T2$mean
        vars_T2  <- preds_T2$var
      } else if (method == "thinnedSV") {
        preds_T2 <- predictions_scaled_thinned(fit = fit,
                                               locs_pred = test_T2_X,
                                               m = m.pred,
                                               nsims = 50)
        means_T2 <- preds_T2$means
        vars_T2  <- apply(preds_T2$samples, 1, var)
      } else if (method == "twingp") {
        preds_T2 <- twingp::twingp(trainX, trainY, test_T2_X)
        means_T2 <- preds_T2$mu
        vars_T2  <- (preds_T2$sigma)^2
      } else if (method == "thinned_twingp_full") {
        Tthin    <- computeThinningNumber(trainX, 20)
        preds_T2 <- thinned_twingp_full(trainX, trainY, test_T2_X, T = Tthin)
        means_T2 <- preds_T2$unweighted$mu
        vars_T2  <- (preds_T2$unweighted$sigma)^2
      } else if (method == "thinned_lagp") {
        Tthin    <- computeThinningNumber(trainX, 20)
        preds_T2 <- thinned_lagp(trainX, trainY, test_T2_X, T = Tthin)
        means_T2 <- preds_T2$mu
        vars_T2  <- (preds_T2$sigma)^2
      } else if (method == "tempgp") {
        fit_tmp <- tempGP(as.matrix(trainX), as.vector(trainY),
                          fast_computation    = FALSE,
                          max_thinning_number = 30L,
                          limit_memory        = 5000L)
        preds   <- predict(fit_tmp, as.matrix(test_T2_X))
        means_T2 <- preds
        vars_T2  <- NULL   
      } else if (method == "SVecchia+t") {
  
        idx_train   <- 1:nrow(trainX)
        idx_test_T2 <- (nrow(trainX) + 1):(nrow(trainX) + nrow(test_T2_X))
        test_T2_X_t <- cbind(idx_test_T2, test_T2_X)
        preds_T2 <- predictions_scaled(fit = fit,
                                       locs_pred = test_T2_X_t,
                                       m = m.pred,
                                       scale = scale,
                                       nsims = 50)
        means_T2 <- preds_T2$means
        vars_T2  <- apply(preds_T2$samples, 1, var)
      }
    })["elapsed"]
    runtime_T2 <- as.numeric(runtime_T2)
    
    ##RMSE, NLPD
    RMSE_T2 <- sqrt(mean((t2$normPW - means_T2)^2))
    
    if (is.null(vars_T2)) {
      NLPD_T2 <- NA_real_   # for tempgp
    } else {
      NLPD_T2 <- 0.5 * mean(((t2$normPW - means_T2)^2) / vars_T2 +
                              log(2 * pi * vars_T2))
    }
  
    results_df <- bind_rows(
      results_df,
      data.frame(
        Dataset     = dataset,
        Method      = method,
        RMSE_T2     = RMSE_T2,
        Runtime_Fit = fit_runtime,
        Runtime_T2  = runtime_T2,
        NLPD_T2     = NLPD_T2,
        stringsAsFactors = FALSE
      )
    )
  }
}

# --------- MAKE TABLE 8 ---------

wt_map <- c(
  "Inland Wind Farm Dataset2(WT1)"  = "WT1",
  "Inland Wind Farm Dataset2(WT2)"  = "WT2",
  "Offshore Wind Farm Dataset2(WT3)"= "WT3",
  "Offshore Wind Farm Dataset2(WT4)"= "WT4"
)
results_df$WT <- wt_map[results_df$Dataset]

method_label_map <- c(
  "H-laGP"              = "laGP",
  "twingp"              = "twinGP",
  "SVecchia"            = "SV(x)",
  "SVecchia+t"          = "SV(x,t)",
  "thinned_lagp"        = "thinned laGP",
  "thinned_twingp_full" = "thinned twinGP",
  "thinnedSV"           = "thinned SV",
  "tempgp"              = "tempGP"
)

rmse_mat <- with(results_df,
                 tapply(RMSE_T2, list(Method, WT), mean))

nlpd_mat <- with(results_df,
                 tapply(NLPD_T2, list(Method, WT), mean))

rmse_mat <- rmse_mat[, c("WT1","WT2","WT3","WT4"), drop = FALSE]
nlpd_mat <- nlpd_mat[, c("WT1","WT2","WT3","WT4"), drop = FALSE]

rmse_avg <- rowMeans(rmse_mat, na.rm = TRUE)
nlpd_avg <- rowMeans(nlpd_mat, na.rm = TRUE)

table8_df <- data.frame(
  Method      = method_label_map[rownames(rmse_mat)],
  RMSE_WT1    = rmse_mat[,"WT1"],
  RMSE_WT2    = rmse_mat[,"WT2"],
  RMSE_WT3    = rmse_mat[,"WT3"],
  RMSE_WT4    = rmse_mat[,"WT4"],
  RMSE_Avg    = rmse_avg,
  NLPD_WT1    = nlpd_mat[,"WT1"],
  NLPD_WT2    = nlpd_mat[,"WT2"],
  NLPD_WT3    = nlpd_mat[,"WT3"],
  NLPD_WT4    = nlpd_mat[,"WT4"],
  NLPD_Avg    = nlpd_avg,
  row.names   = NULL
)

table8_disp <- table8_df
table8_disp[,-1] <- round(table8_disp[,-1], 2)



# --------- MAKE TABLE 9 ---------

results_df$Runtime_Total_min <- (results_df$Runtime_Fit +
                                   results_df$Runtime_T2) / 60

runtime_by_method <- aggregate(Runtime_Total_min ~ Method,
                               data = results_df, FUN = mean)


all_method_labels <- method_label_map[
  c("H-laGP","twingp","SVecchia","SVecchia+t",
    "thinned_lagp","thinned_twingp_full","thinnedSV","tempgp")
]

runtime_row <- setNames(
  as.list(rep(NA_real_, length(all_method_labels))),
  all_method_labels
)

for (i in seq_len(nrow(runtime_by_method))) {
  m  <- runtime_by_method$Method[i]
  ml <- method_label_map[m]
  runtime_row[[ml]] <- runtime_by_method$Runtime_Total_min[i]
}

table9_df <- data.frame(
  Dataset6 = "Dataset 6",
  t(as.data.frame(runtime_row)),
  row.names = NULL
)


table9_disp <- table9_df
table9_disp[,-1] <- round(as.numeric(unlist(table9_disp[,-1])), 1)


table8_disp
table9_disp

#write.csv(table8_disp, "Table8_DSWE6.csv", row.names = FALSE)
#write.csv(table9_disp, "Table9_Runtime_DSWE6.csv", row.names = FALSE)
