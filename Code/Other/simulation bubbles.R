library(exuber)
library(psymonitor)
library(rugarch)
library(devtools)
library(rbmi)
library(dplyr)
library(dqshiny)
library(ggplot2)
library(ghyp)
library(lubridate)
library(strucchange)
library(changepoint)
library(Rssa)
library(tidyr)
library(qcc)
library(cpm)
library(Rssa)
library(lattice)
library(fANCOVA)
library(stats)
library(HoRM)
library(WaveletComp)
#install_url('https://cran.r-project.org/src/contrib/Archive/psymonitor/psymonitor_0.0.2.tar.gz')
#install_url('https://cran.r-project.org/src/contrib/Archive/dqshiny/dqshiny_0.0.4.tar.gz')
add_class <- function(x, ...) {
  class(x) <- append(c(...), class(x))
  x
}
#DATA
#https://rdrr.io/rforge/rugarch/man/sp500ret.html
#https://rdrr.io/rforge/rugarch/man/dmbp.html
#data(dmbp)
data(sp500ret)
#auto.arima(sp500ret$SP500RET)
data_garch = sp500ret*100
#GARCH
spec <- ugarchspec(variance.model = list(model = "eGARCH",garchOrder = c(1,1))
                   ,distribution.model="ghyp")
fit = ugarchfit(data = data_garch , spec = spec)

#residuals
garch_resid <- residuals(fit, standardize = FALSE)
# extract standardized residuals
garch_resid_std <- residuals(fit, standardize = TRUE)
# standardized residuals to be used in ugarchsim()
custom_dist = list(
  name = "sample",
  distfit = matrix(garch_resid_std,  ncol = 1)
)

m_ <- fit@model$maxOrder

garch_sim <- ugarchsim(fit, 
                       n.sim = length(garch_resid_std),
                       m.sim = 1, 
                       presigma = tail(fit@fit$sigma, m_),
                       prereturns = tail(sp500ret, m_),
                       # preresiduals ARE NOT standardized:
                       preresiduals = tail(garch_resid, m_),
                       startMethod = "sample",
                       # distfit in custom.dist ARE standardized
                       custom.dist = custom_dist)

# extract simulated series
sp_sim <- garch_sim@simulation$seriesSim
# test
sp_df <- cbind(sp500ret, sim = sp_sim) %>% 
  setNames(c("ret", "sim")) %>% 
  as_tibble(rownames = "date") %>% 
  mutate(date = as.Date(date),
         diff = sim - ret)

egarch<-function(n_sim) {
  #custom_dist = list(
  # name = "sample",
  #distfit = matrix(garch_resid_std[1:n_sim,],  ncol = 1)
  #)
  
  garch_sim <- ugarchsim(fit, 
                         n.sim = n_sim,
                         m.sim = 1, 
                         presigma = tail(fit@fit$sigma, m_),
                         prereturns = tail(sp500ret, m_),
                         # preresiduals ARE NOT standardized:
                         preresiduals = tail(garch_resid, m_),
                         startMethod = "sample",
                         # distfit in custom.dist ARE standardized
                         #custom.dist = custom_dist
                         )
  values = as.numeric(garch_sim@simulation$seriesSim)
  #remove(custom_dist)
  return(values)
}


bubble <- function(n, te = 0.4 * n, tf = te + 0.2 * n , tr = tf + 0.1*n,
                   c = 1, c1 = 1, c2 = 1, eta = 0.6, alpha = 0.6, beta = 0.5) {
  
  drift <- c*n^(-eta)
  delta <- 1 + c1 * n^(-alpha)
  gamma <- 1 - c2 * n^(-beta)
  y <- 100
  err <- egarch(n)
  
  for (t in 2:n) {
    if (t < te) {
      y[t] <- drift + y[t - 1] + err[t]
    } else if (t >= te & t <= tf) {
      y[t] <- delta * y[t - 1] + err[t]
    } else if (t > tf & t <= tr ) {
      y[t] <- gamma * y[t - 1] + err[t]
    } else {
      y[t] <- drift + y[t - 1] + err[t]
    }
  }
  y %>%
    add_class("sim") 
}


bubble_2 <- function(n,
                    te1 = 0.3 * n, tf1 = te1 + 0.1 * n , tr1 = tf1 + 0.1*n,
                    te2 = 0.8 * n, tf2 = te2 + 0.1 * n , tr2 = tf2 + 0.1*n,
                    c = 1, c1 = 1, c2 = 1, eta = 0.6, alpha = 0.6, beta = 0.7) {
  
  drift <- c*n^(-eta)
  delta <- 1 + c1 * n^(-alpha)
  gamma <- 1 - c2 * n^(-beta)
  y <- 100
  err <- egarch(n)
  
  for (t in 2:n) {
    if (t < te1) {
      y[t] <- drift + y[t - 1] + err[t] # normal
    } else if (t >= te1 & t <= tf1) {
      y[t] <- delta * y[t - 1] + err[t] # bubble1
    } else if (t > tf1 & t <= tr1 ) {
      y[t] <- gamma * y[t - 1] + err[t] # collapse 1
    }  else if (t > tr1 + 1 & t < te2) {
      y[t] <- drift + y[t - 1] + err[t] # normal 2
    }  else if (t >= te2 + 1 & t <= tf2) {
      y[t] <- delta * y[t - 1] + err[t] # bubble 2
    }  else if (t > tf2 + 1 & t <= tr2) {
      y[t] <- gamma * y[t - 1] + err[t] # collapse 2
    } else {
      y[t] <- drift + y[t - 1] + err[t] # normal 3
    }
  }
  y %>%
    add_class("sim")
}

set.seed(000)
time = 100
disturbing <- bubble(time, te = time*0.4, tf= time*0.6, tr = time*0.7, beta = 0.5, alpha = 0.6)
sudden <- bubble(time, time*0.3, tf= time*0.5, tr = time*0.51, beta = 0.15, alpha = 0.6, eta = 0.03)
smooth <- bubble(time, time*0.4, tf= time*0.6, tr = time*0.8, beta = 0.9, alpha = 0.6)
double<-bubble_2(time, alpha = 0.6, beta = 0.7, eta = 0.4)

autoplot(disturbing)
autoplot(smooth)
autoplot(sudden)
autoplot(double)

plot(diff(log(smooth)),type='l')
plot(diff(log(disturbing)),type='l')
plot(diff(log(sudden)),type='l')
plot(diff(log(double)),type='l')

plot(abs(diff(log(smooth))),type='l')
plot(abs(diff(log(disturbing))),type='l')
plot(abs(diff(log(sudden))),type='l')
plot(abs(diff(log(double))),type='l')


#TEST
sp500ret$data <- row.names(sp500ret)
sp500ret$data <- ymd(sp500ret$data)

series <- abs(diff(log(disturbing)))
              
test <- processStream(series,"Kolmogorov-Smirnov",ARL0=10000, startup=20)
plot(series, type = "l", xlab = "Observation", ylab = "", bty = "l")
abline(v = test$detectionTimes, col = 'blue')
abline(v = test$changePoints, lty = 2, col = 'red')
test

ewma <- numeric(length(series))
ewma[1] <- series[1]^2
lambda <- 0.94
for (i in 2:length(ewma)){
  ewma[i] <- lambda * ewma[i - 1] + (1 - lambda) * series[i]^2
}
plot(ewma, type = "l", xlab = "Observation", ylab = "", bty = "l")

abline(v = test$changePoints, lty = 2)

radf_sim <- radf(disturbing)
autoplot(radf_sim)
datestamp(radf_sim)

#REGRESSOGRAM, LOCAL POLYNOMIAL AND WAVELET
#reg<-regressogram(1:length(series),series, nbins = 254)
#plot(reg)
l<-loess.as(1:length(series),series, degree = 2, criterion = 'aicc', plot = TRUE, family = 'symmetric')
my.w <- analyze.wavelet(as.data.frame(series),'series',
                        loess.span = 0,
                        dt = 1,
                        lowerPeriod = 4,
                        upperPeriod = 800,
                        #dj = 1/10,
                        make.pval = TRUE, n.sim = 10)

wt.image(my.w, color.key = "quantile", n.levels = 100,
         legend.params = list(lab = "wavelet power levels", mar = 4.7))

rrec<-reconstruct(my.w, plot.waves = FALSE, lwd = c(1,2),
            legend.coords = "topleft")

res_w<-(series-rrec$series$series.r)
plot(res_w, type='l')

#chart<-ewma(rrec$series$series.r, lambda=0.94, center = mean(abs(res_w)), std.dev =  sd(abs(res_w)))
#chart$violations

#SIMULATIONS
matrix_empty = matrix(nrow = 100, ncol = 100)

set.seed(002)
cyc = c(1:100)
DIS <- data.frame(matrix_empty)
for (i in cyc){
  DIS[,i] <- bubble(time, te = time*0.4, tf= time*0.6, tr = time*0.7, beta = 0.5, alpha = 0.6)
}

results_DIS <- list()
ADF_DIS <- list()
for (i in cyc){
  ss <- abs(diff(log(DIS[,i])))
  test <- processStream(ss,"Kolmogorov-Smirnov",ARL0=10000, startup=20)
  results_DIS[[i]]<-test$changePoint
  remove(test)
  radf_sim <- radf(DIS[,i])
  d<-datestamp(radf_sim)
  ADF_DIS[[i]] <- c(d$series1$Start,d$series1$End)
  remove(d)
  remove(radf_sim)
  print(i)
}
results_DIS
ADF_DIS
correct<-results_DIS[lengths(results_DIS) == 3]


set.seed(003)
SUD <- data.frame(matrix_empty)
for (i in cyc){
  SUD[,i] <- c(bubble(time, time*0.3, tf= time*0.5, tr = time*0.51, beta = 0.15, alpha = 0.6))
}

results_SUD <- list()
ADF_SUD <- list()
for (i in cyc){
  ss <- abs(diff(log(SUD[,i])))
  test <- processStream(ss,"Kolmogorov-Smirnov",ARL0=10000, startup=20)
  results_SUD[[i]]<-test$changePoints
  remove(test)
  radf_sim <- radf(SUD[,i])
  d<-datestamp(radf_sim)
  ADF_SUD[[i]] <- c(d$series1$Start,d$series1$End)
  remove(d)
  remove(radf_sim)
  print(i)
}
results_SUD
correct<-results_SUD[lengths(results_SUD) == 2]
ADF_SUD

set.seed(004)
SMO <- data.frame(matrix_empty)
for (i in cyc){
  SMO[,i] <- c(bubble(time, time*0.4, tf= time*0.6, tr = time*0.8, beta = 0.9, alpha = 0.6))
}

results_SMO <- list()
ADF_SMO <- list()
for (i in cyc){
  ss <- abs(diff(log(SMO[,i])))
  test <- processStream(ss,"Kolmogorov-Smirnov",ARL0=10000, startup=20)
  results_SMO[[i]]<-test$changePoints
  remove(test)
  radf_sim <- radf(SMO[,i])
  d<-datestamp(radf_sim)
  ADF_SMO[[i]] <- c(d$series1$Start,d$series1$End)
  remove(d)
  remove(radf_sim)
  print(i)
}
results_SMO
correct<-results_SMO[lengths(results_SMO) == 3]
ADF_SMO

#set.seed(005)
#DUB <- data.frame(matrix_empty)
#for (i in cyc){
#  DUB[,i] <- c(bubble_2(time, alpha = 0.7, beta = 0.8, eta = 0.4))
#}

#results_DUB <- list()
#for (i in cyc){
#  ss <- abs(diff(log(DUB[,i])))
#  test <- processStream(ss,"Kolmogorov-Smirnov",ARL0=10000, startup=20)
#  results_DUB[[i]]<-test$changePoints
#  remove(test)
#}
#results_DUB

#S&P500
library(bbdetection)
SP500 <- as.data.frame(sp500m)
SP500$Date <- row.names(SP500)
SP500$Date <- ymd(SP500$Date)
ggplot(SP500, aes(Date, SP500))+geom_line()
series <- abs(diff(log(SP500$SP500)))
plot(SP500$Date[2:840],series, type='l')

test <- processStream(series,"Kolmogorov-Smirnov", ARL0=2000, startup=24)
breack<-SP500$Date[test$changePoints]
plot(SP500$Date, SP500$SP500, type = "l", xlab = "Observation", ylab = "", bty = "l")
#abline(v = test$detectionTimes, col = 'blue')
abline(v = breack, lty = 2, col = 'red')
test

ewma <- numeric(length(series))
ewma[1] <- series[1]^2
lambda <- 0.94
for (i in 2:length(ewma)){
  ewma[i] <- lambda * ewma[i - 1] + (1 - lambda) * series[i]^2
}
plot(SP500$Date[2:840], ewma, type = "l", xlab = "Observation", ylab = "", bty = "l")
abline(v = breack, lty = 2, col = 'red')

#s <- ssa(SP500$SP500, L=2)
#plot(s)
#plot(s, 'series')
#plot(s, 'vectors')
#plot(wcor(s))
#wcor(s)

#gr <- grouping.auto(s, grouping.method = "wcor", nclust=2)
#rec <- Rssa::reconstruct(s, groups = gr)

#specs <-lapply(rec, function(x) spectrum(x, plot = FALSE)$spec)
#w.tree <- seq(0, length.out = length(rec$'1'),
#             by = length(SP500$SP500)) 
#xyplot(rec$'1'+rec$'2'~ w.tree, data = specs,
#       superpose = FALSE, type = "l", xlab = NULL, ylab = NULL,
#       auto.key = list(lines = TRUE, points = TRUE,
#                      column = 2))

#plot(rec$'2',type='l')
#plot(SP500$Date, res, type='l')

l<-loess.as(1:length(SP500$SP500),SP500$SP500, degree = 2, criterion = 'aicc', plot = TRUE, family = 'symmetric')
plot(SP500$Date,abs(l$residuals), type='l')

s_2 <- ssa(series)
sequence = c(1:100)
sums <- c(rep(0,length(series)))
HZ <- list(0.038)
for (i in sequence){
  a <- runif(1, 0.1, 0.9)
  group <- grouping.auto(s_2, base = "series",
                         freq.bins = HZ,
                         threshold = a
  )
  result <- Rssa::reconstruct(s_2, group)
  sums <- c(Reduce(`+`, result)+sums)
  print(i)
}
sums<-sums/length(sequence)
plot(SP500$Date[2:840], sums, type='l')
abline(v = breack, lty = 2, col = 'red')

#res<-abs(rec$'2')
#plot(res, type='l')
#s_res<-ssa(res)
#plot(s_res, type = "paired", plot.contrib = FALSE)


#MANUAL
#plot (s, 'series', plot.contrib = FALSE)
s_2 <- ssa(SP500$SP500, L=length(SP500$SP500)/2)
plot(s_2, 'series')
lst <- grouping.auto(s_2, grouping.method = "wcor",
                     nclust=2)
print(lst)
plot(lst)
# Check separability
w <- wcor(s_2, groups = lst)
plot(w)
g <- Rssa::reconstruct(s_2, groups = lst)
plot(g, add.residuals = FALSE, idx = 1:length(lst), plot.method = "xyplot", superpose = FALSE)

SP500$smooth<-g$'1'
SP500$res<-c(0,ewma)*50000-500

ggplot(SP500,aes(Date))+
  geom_line(aes(y=SP500), colour="black")+ geom_line(aes(y=smooth), colour="blue")+geom_line(aes(y=res), color="green") + 
  geom_vline(xintercept=breack, linetype='dashed', color='red') + 
  labs(title = 'S&P 500 speculative bubble') + scale_color_identity(name= NULL, 
                                                                 labels = c(black = "S&P500", green = "EWMA", blue = "SSA"), guide='legend')
  
crit_values<-radf_mc_cv(nrow(SP500), nrep = 2000, seed=123)
results<-radf(SP500$SP500)
datestamp(results, crit_values)
c95<-as.data.frame(crit_values[["bsadf_cv"]])

a<-autoplot(results, crit_values) + labs(title = "S&P 500")

#Bayesian change point
library(Rbeast)
s_3 <- ssa(SP500$SP500, L=length(SP500$SP500)/2)
plot(s_3, 'series')
lst <- grouping.auto(s_2, grouping.method = "wcor",
                     nclust=50)
print(lst)
plot(lst)
# Check separability
w <- wcor(s_3, groups = lst)
plot(w)
g <- Rssa::reconstruct(s_3, groups = lst)
plot(g, add.residuals = FALSE, idx = 1:length(lst), plot.method = "xyplot", superpose = FALSE)

SP500$smooth_2<-g$'1'

out = beast(series, season = 'none', precPriorType = 'uniform', mcmc.chains = 10, mcmc.samples = 10000, mcmc.seed = 001, mcmc.burning  = 1000)
plot(out)
change_p<-SP500$Date[out[["trend"]][["cp"]]]
ggplot(SP500,aes(Date))+
  geom_line(aes(y=SP500), colour="black")+ geom_line(aes(y=smooth_2), colour="blue")+geom_line(aes(y=res), color="green") + 
  geom_vline(xintercept=change_p, linetype='dashed', color='red') + 
  labs(title = 'S&P 500 speculative bubble') + scale_color_identity(name= NULL, 
                                                                    labels = c(black = "S&P500", green = "EWMA", blue = "SSA"), guide='legend')

series<-abs(diff(log(sudden)))
out = beast(series, season = 'none', precPriorType = 'uniform', mcmc.chains = 10, mcmc.samples = 10000, mcmc.seed = 002, mcmc.burning  = 1000)
plot(out)

series<-abs(diff(log(disturbing)))
out = beast(series, season = 'none', precPriorType = 'uniform', mcmc.chains = 10, mcmc.samples = 10000, mcmc.seed = 003, mcmc.burning  = 1000)
plot(out)

series<-abs(diff(log(smooth)))
out = beast(series, season = 'none', precPriorType = 'uniform', mcmc.chains = 10, mcmc.samples = 10000, mcmc.seed = 003, mcmc.burning  = 1000)
plot(out)

############
#SSA
s <- ssa(series)
plot(s)
plot (s, 'series', plot.contrib = FALSE)
plot (s, 'vectors')
plot(wcor(s))
plot(s, type = "paired")

#AUTOMATIC SELECTION
sequence = c(1:100)
sums <- c(rep(0,length(series)))
HZ <- list(0.033)
for (i in sequence){
  a <- runif(1, 0.1, 0.9)
  group <- grouping.auto(s, base = "series",
                         freq.bins = HZ,
                         threshold = a
  )
  result <- Rssa::reconstruct(s, group)
  sums <- c(Reduce(`+`, result)+sums)
  print(i)
}
sums<-sums/length(sequence)
plot(sums, type='l')
abline(h=mean(sums), col="blue", lwd=1, lty=1)


#MANUAL
#plot (s, 'series', plot.contrib = FALSE)
lst <- grouping.auto(s, grouping.method = "wcor",
                     nclust=5)
print(lst)
plot(lst)
# Check separability
w <- wcor(s, groups = lst)
plot(w)
g <- Rssa::reconstruct(s, groups = lst)
plot(g, add.residuals = FALSE, idx = 1:length(lst), plot.method = "xyplot", superpose = FALSE)


#plot (s, 'series', plot.contrib = FALSE)
#g <- grouping.auto(s, base = "series",
#                    freq.bins = HZ,
#                   threshold = 0.9
#                    )

#plot(g, order = TRUE, type = "b", idx = 1:10)
#r1 <- reconstruct(s, g)
#plot(r1, plot.method = "xyplot", superpose = TRUE, add.residuals = FALSE)
#plot(r1, add.residuals = FALSE)
#plot(Reduce(`+`, r1), type='l')

specs <-lapply(g, function(x) spectrum(x, plot = FALSE)$spec)
w.tree <- seq(0, length.out = length(g$'1'),
              by = length(series)/100) 
xyplot(g$'1'+g$'2'+g$'3'+g$'4'~ w.tree, data = specs,
       superpose = FALSE, type = "l", xlab = NULL, ylab = NULL,
       auto.key = list(lines = TRUE, points = FALSE,
                       column = 4))

model <- sums
plot(model, type = 'l')
res <- (series-sums)
#res <- series - c(sums/length(sequence))
plot(res, type='l')
mean(abs(res))

confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

#xyplot(series + sums ~ sp500ret$data, type = "l") 


spec.pgram(SP500$SP500, detrend = FALSE)
#spec.pgram(series, detrend = FALSE, log='no')
#test <- cpt.mean(sums, method = "PELT")
#plot(test, type = "l", cpt.col = "blue", xlab = "Index", cpt.width = 4)
#cpts(test)
#dat <- tibble(ylag0 = abs(res),
#              ylag1 = lag(abs(res))

#) %>%
#  drop_na()


#cusum <- Fstats(ylag0 ~ ylag1, data = dat)
#plot(cusum)

## SD-EWMA
#data_ew <- data.frame(timestamp = 1:length(series), value = series)

#result <- OcpSdEwma(
#  data = sums,
#  n.train = 10,
#  threshold = 0.01,
#  l = 3
#)
#plot <- cbind(data_ew, result)
#PlotDetections(plot, title = "SD-EWMA ANOMALY DETECTOR")

plot(series, type='l', col='blue', lty = 1,pch=19,lwd=1)
lines(sums, col='green', lty = 1,pch=19,lwd=1)


#SIMULATION TEST

train <- function(n, c = 1, eta = 0.6, drift = 1) {
  y <- 100
  err <- egarch(n)
  
  for (t in 2:n) {
    y[t] <- eta + drift*y[t - 1] + err[t]
  }
  y %>%
    add_class("sim") 
}

set.seed(987)
sequence_2 = c(1:100)
sums_2 <- c(rep(0,length(series)))
max <- c(0)
for (i in sequence_2){
  data_train <- train(length(series))
  sums_2 <- c(data_train + sums_2)
  print(i)
  max <- max(abs(diff(log(data_train)))) + max
}

sums_2 <- sums_2/length(sequence_2)
plot(abs(diff(log(sums_2))), type='l')
abline(h=mean(sums_2), col="blue", lwd=1, lty=1)
max<-max/length(sequence_2)
max



######
s_1<-grouping.auto.wcor(s, groups = 1:5, nclust=4)
plot(s_1)

recon <- reconstruct(s, groups = list(1))
plot(recon)

res1 <- reconstruct(s, groups = list(1))
trend = res1$F1
plot(trend,type='l')

res.trend<- residuals(res1)
spec.pgram(res.trend, detrend = FALSE, log = "no")
s2 <- ssa(res.trend, L=50)
plot(s2)
plot(s2, type = "paired", idx = 1:9, plot.contrib = FALSE)
#res2 <- reconstruct(s2, groups=list(1:2))
#seasonality <- res2$F1
#res <- residuals(res2)
#plot(res1, add.residuals = FALSE, col = c("black", "red"))
plot(s2, 'series')
plot(wcor(s2))

res2 <- reconstruct(s2, groups = list(1:4))
volatility <- res2$F1

res.volatility <- residuals(res2)
spec.pgram(res.volatility, detrend = FALSE, log = "no")
s.env <- ssa(res.volatility^2, L=30)
rsd <- sqrt(reconstruct(s.env, groups=list(1))$F1)

#plot(res, type='l')
#lines(rsd, type='l', col = "red")
#lines(-rsd, type='l', col = "red")

plot(volatility, type ='l')
plot(trend, type ='l')



