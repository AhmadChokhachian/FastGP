#Please switch "path" with the actual download direction
setwd("path/Dataset6")

datasets <- c("Offshore Wind Farm Dataset2(WT3)")

results_df <- data.frame(Dataset = character(), Method = character(),
                         Seed = integer(), RMSE_T2 = numeric())

m.pred = 140
m.hybrid = 30
m.est = 30

methods <- c('thinnedSV','Svecchia','twingp','thinned_twingp_full')


for (dataset in datasets) {
  df <- read.csv(paste0(dataset, ".csv"), header = TRUE, sep = ",")
  
  df$time <- ymd_hms(df$time)
  
  df$wind_direction_sin <- sin(df$D)
  df$wind_direction_cos <- cos(df$D)
  
  xCol = c(2, 3, 4, 5, 6, 7, 9, 10) 
  yCol = c(8) 
  
  if (dataset %in% c("Offshore Wind Farm Dataset2(WT3)", "Offshore Wind Farm Dataset2(WT4)")) {
    t1 <- subset(df, year(time) %in% c(2007, 2008))
    t2 <- subset(df, year(time) == 2009)
  } else {
    t1 <- subset(df, year(time) %in% c(2008, 2009))
    t2 <- subset(df, year(time) == 2010)
  }

  trainX <- as.matrix(t1[, xCol])
  trainY <- as.numeric(t1[, yCol])
  test_T2_X <- as.matrix(t2[, xCol])

  for (seed in 1:10) {
    set.seed(seed)  
    
    for (method in methods) {
      print(paste("Dataset:", dataset, "Method:", method, "Seed:", seed))
      
      fit_runtime <- system.time({
        if (method == 'Svecchia') {
          scale = 'parms'
          fit <- fit_scaled(trainY, trainX, ms = m.est, scale = scale,
                            trend = 'linear', print.level = 2, find.vcf = FALSE)
        } else if (method == 'thinnedSV') {
          T <- computeThinningNumber(trainX, 20)
          fit <- fit_scaled_thinned(trainY, trainX, ms = m.est,trend = 'zero',scale = 'parms',max.it = 64, 
                                    thinnedBins = createThinnedBins(trainX, trainY, T), T = T)
        } else if (method == 'twingp') {
          fit <- list(trainX = trainX, trainY = trainY) 
        } else if (method == 'thinned_twingp_full') {
          fit <- list(trainX = trainX, trainY = trainY) 
        }
      })['elapsed']
      
      runtime_T2 <- system.time({
        if (method == 'Svecchia') {
          preds_T2 <- predictions_scaled(fit = fit, locs_pred = test_T2_X, 
                                         m = m.pred, scale = scale, nsims = 50)
          means_T2 <- preds_T2$means
        } else if (method == 'thinnedSV') {
          preds_T2 <- predictions_scaled_thinned(fit = fit, locs_pred = test_T2_X, 
                                                 m = m.pred)
          means_T2 <- preds_T2
        } else if (method == 'twingp') {
          preds_T2 <- twingp::twingp(trainX, trainY, test_T2_X)
          means_T2 <- preds_T2$mu
        } else if (method == 'thinned_twingp_full') {
          T <- computeThinningNumber(trainX, 20)
          preds_T2 <- thinned_twingp_full(trainX, trainY, test_T2_X, T = T)
          means_T2 <- preds_T2$unweighted$mu
        } 
      })['elapsed']
      
      RMSE_T2 <- sqrt(mean((t2$normPW - means_T2)^2))

      results_df <- bind_rows(results_df, data.frame(
        Dataset = dataset, Method = method, Seed = seed, RMSE_T2 = RMSE_T2
      ))
      
      
     }
  }
}

library(ggplot2)
results_df$Method <- factor(
  results_df$Method,
  levels = c("twingp", "thinned_twingp_full", "Svecchia", "thinnedSV"),
  labels = c("twinGP", "thinned twinGP", "SV(x)", "thinnedSV")
)
cols <- c(
  "twinGP"         = "#F7EA5D",
  "thinned twinGP" = "#273B78",
  "SV(x)"          = "#B569FF",
  "thinnedSV"      = "#BFBFBF"
)
stability_plot <- ggplot(results_df, aes(x = Method, y = RMSE_T2, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = cols) +
  labs(x = "", y = "RMSE") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

stability_plot
