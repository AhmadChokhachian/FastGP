#Please switch "path" with the actual download direction
setwd("path/robot_data")

plot_custom_pacf <- function(file_LHS, file_mod, file_high, lags = 35) {

  data_LHS <- read.csv(file_LHS)[[1]]
  data_mod <- read.csv(file_mod)[[1]]
  data_high <- read.csv(file_high)[[1]]

  pacf_LHS <- pacf(data_LHS, lag.max = lags, plot = FALSE)
  pacf_mod <- pacf(data_mod, lag.max = lags, plot = FALSE)
  pacf_high <- pacf(data_high, lag.max = lags, plot = FALSE)
  
  max_y <- max(
    abs(c(pacf_LHS$acf, pacf_mod$acf, pacf_high$acf))
  )
  y_limit <- c(-max_y, max_y)
  
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))  
  pacf_plot <- function(data, title_text, ylab = FALSE) {
    pacf(data, lag.max = lags, main = title_text,
         ylim = y_limit,
         ylab = if (ylab) "Partial ACF" else "",
         cex.main = 2, cex.lab = 1.8, cex.axis = 1.5)
  }
  

  pacf_plot(data_LHS, "T=2", ylab = TRUE)
  pacf_plot(data_mod, "T=15", ylab = FALSE)
  pacf_plot(data_high, "T=23", ylab = FALSE)
}


plot_custom_pacf(
  file_LHS = "low_autocorrelation/rep_1_train.csv",
  file_mod = "medium_autocorrelation/rep_1_train.csv",
  file_high = "high_autocorrelation/rep_1_train.csv"
)



