#-----------------------------------------------------------
# Path
root_path = "/Users/louissellier/Desktop/R_directory1/R_models/MA_model"
data_path = file.path(root_path, "data")
data_file = "data.csv"


#-----------------------------------------------------------
# Data
df = read.csv(file.path(data_path, data_file))
colnames(df) <- c('date', 'p')

MA_orders = c(1, 5, 10, 15)

build_ma = function(x, n) {
  out =  matrix(0, length(x), 1)
  for (i in n:length(x)) {
    out[i, 1] = mean(x[(i - n + 1):i])
  }
  colnames(out) = paste0("MA(", as.character(n), ")")
  return(out)
}

build_mat_ma = function(y, n) {
  out = do.call(cbind , lapply(X = n, function(x) build_ma(y, x)))
  return(out)
}

x = build_mat_ma(df$p, MA_orders)[MA_orders[length(MA_orders)]:nrow(df), ]
df = df[MA_orders[length(MA_orders)]:nrow(df), ]

build_ma_diffs = function(x, type = 1) {
  out = NULL
  names = NULL
  if (type == 1) {
    for (i in 1:(ncol(x) - 1)) {
      out = cbind(out, x[ , i] - x[ , i + 1])
      names = c(names, paste0(colnames(x)[i]," - ", colnames(x)[i + 1]))
    }
  }
  else if (type == 2) {
    for (i in 2:ncol(x)) {
      out = cbind(out, x[ , 1] - x[ , i])
      names = c(names, paste0(colnames(x)[1]," - ", colnames(x)[i]))
    }
  }
  colnames(out) = names
  return(out)
}

x = build_ma_diffs(x, type = 2)
df = data.frame(df, x)
p_diff = cbind(diff(df$p))
colnames(p_diff) = "p_diff"
df = data.frame(df[2:nrow(df), ], p_diff)


#-----------------------------------------------------------
# Calibrate model
in_sample = 800
dfin = df[1:in_sample, ]

Y = cbind(dfin$p_diff)
X = as.matrix(dfin[ ,3:(ncol(df) - 1)])

build_ols_coefficients = function(X,Y) {
  return(solve(t(X) %*% X) %*% t(X) %*% Y)
}

coef = build_ols_coefficients(X, Y)

Yf = X %*% coef
E = Yf - Y
dfin = data.frame(dfin, Yf, E)
View(dfin)

Rsqr = (sd(Y) - sd(E)) / sd(Y) 
Updown = Y/abs(Y) - Yf/abs(Yf)
accuracy = sum(as.numeric(Updown == 0), na.rm = T)/nrow(Updown)

plot_generator = function(Y, Yf, E, x_axis = NULL) {
  
  # create a default x axis
  if (is.null(x_axis)) {
    x_axis = matrix(1:nrow(Y),nrow(Y),1)
  }
  # create a dataframe
  DF = data.frame(x_axis,Y,Yf,E)
  # load plotting libraries
  library(ggplot2)
  library(gridExtra)
  
  # plot forecasted versus observed timeseries
  plot_1 = ggplot(data = DF, aes(x = x_axis)) +
    geom_line(aes(y = Y, color = 'Y')) +
    geom_point(aes(y = Y, color = 'Y'), size = 0.5) +
    geom_line(aes(y = Yf, color = 'Yf')) +
    geom_point(aes(y = Yf, color = 'Yf'), size = 0.5) +
    scale_color_manual(values = c('Y' = 'red', 'Yf' = 'blue')) +
    #scale_x_date(date_labels = '%Y',date_breaks = '2 year') +
    labs(color = '') +
    theme(legend.position = 'top') +
    xlab('time_index') + ylab('') 
  grid.arrange(plot_1)
  
  #plot residuals timeseries
  plot_2 = ggplot(data = DF, aes(x = x_axis)) +
    geom_line(aes(y = E, color = 'E')) +
    geom_point(aes(y = E, color = 'E'), size = 0.5) +
    scale_color_manual(values = c('E' = 'black')) +
    #scale_x_date(date_labels = '%Y',date_breaks = '2 year') +
    labs(color = '') +
    theme(legend.position = 'top') +
    xlab('date_index') + ylab('') 
  grid.arrange(plot_2)
  
  #plot residuals histogram
  plot_3 = ggplot(data = DF, aes(x = E)) +
    geom_histogram(bins = 40, fill = "skyblue", color = "blue") +
    xlab('error') + ylab('frequency')
  grid.arrange(plot_3)
  
  #residuals box plot
  plot_4 = ggplot(data = DF, aes(y = E)) +
    geom_boxplot()
  grid.arrange(plot_4)
  
  #plot residuals heteroscedascity
  plot_5 = ggplot(data = DF, aes(x = Y, y = E)) +
    geom_point()
  grid.arrange(plot_5)
}

plot_generator(Y, Yf, E)



#-----------------------------------------------------------
# Backtest OOS
out_sample = 200
dfout = df[(in_sample + 1):min(in_sample + 1 + out_sample, nrow(df)), ]
Y = cbind(dfout$p_diff)
X = as.matrix(dfout[ ,3:(ncol(dfout) - 1)])
Yf = X %*% coef
E = Yf - Y
dfout = data.frame(dfout, Yf, E)
View(dfout)
plot_generator(Y, Yf, E)

Rsqr = (sd(Y) - sd(E)) / sd(Y) 
Updown = Y/abs(Y) - Yf/abs(Yf)
accuracy = sum(as.numeric(Updown == 0), na.rm = T)/nrow(Updown)


#-----------------------------------------------------------
# Backtest In-Out strategy
base = 100
ret = 1 + dfout$p_diff/dfout$p
port = c(base)
sp500 = c(base)
for (i in 1:nrow(dfout)) {
  if (dfout$Yf[i] < 0) {
    port = rbind(port, port[i - 1])
  }
  else {
    port = rbind(port, port[i - 1] * ret[i])
  }
  sp500 = rbind(sp500, sp500[i - 1] * ret[i])
}

df_strat = data.frame(port, sp500, x_axis = 1:nrow(port)) 
ggplot(data = df_strat, aes(x = df_strat$x_axis)) + 
  geom_line(aes(y = df_strat$port), color = "blue") +
  geom_line(aes(y = df_strat$sp500), color = "red") 


return = port[nrow(port)]/base - 1
gmr = (port[nrow(port)]/base)^(1/nrow(port)) - 1
