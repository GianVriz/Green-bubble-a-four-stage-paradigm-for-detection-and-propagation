#LIBRARIES AND DATA
renixx <- read.csv("~/GitHub/PhD-Thesis/Poject paper 1/renixx.csv")
library(ggplot2)
library(exuber)
library(MultipleBubbles)
#library(tidyverse)
library(dplyr)
library(lubridate)
library(tidyr)
library(vars)
library(astsa)
library(lme4)
library(stargazer)
library(sjPlot)
library(quantreg)
library(nnet)

#DESCRIPTIVE STATISTICS
renixx$date<-ymd(renixx$date)
ggplot(renixx, aes(date, close))+geom_line()

df <- renixx %>% dplyr::select(1,3)
df$log <- log(df$close)
#df_2<-df%>%slice(1:600)
#seed(123)
#crit_values <- radf_mc_cv(nrow(df), seed = 001)
#results<-radf(df$log)
#autoplot(results, crit_values)+labs(title = "Renixx BSADF")


#MONTH GSADF
renixx_month <- renixx %>% mutate(month = format(date, "%m"), year = format(date, "%Y")) %>% group_by(year,month) %>% summarize(Mean_m = mean(close)) %>% as.data.frame
renixx_m <- renixx_month %>% mutate('date' = make_date(year = year, month = month))
remove(renixx_month)
renixx_m$logstock <- log(renixx_m$Mean_m)
ggplot(renixx_m, aes(date, Mean_m))+geom_line()

df_2 <- renixx_m %>% dplyr::select(4,5)
crit_values<-radf_mc_cv(nrow(df_2), nrep = 2000, seed=123)
results<-radf(df_2)
datestamp(results, crit_values)
c95<-as.data.frame(crit_values[["bsadf_cv"]])

autoplot(results, crit_values) + labs(title = "Renixx GSADF")

autoplot(results,crit_values) + 
  labs(title = "GSADF and Renixx index", xlab= "Date", ylab="GSADF test") +
  geom_line(aes(x = date, y = Mean_m/50),data = renixx_m, color = 'black', inherit.aes = FALSE) +
  scale_y_continuous(name = "Renixx index", sec.axis = sec_axis(~.*50, name="Second Axis"))

autoplot(results) + 
  labs(title = "GSADF and Renixx index", xlab= "Date", ylab="GSADF test") +
  geom_line(aes(x = date, y = logstock*1.7),data = renixx_m, color = 'black', inherit.aes = FALSE) +
  scale_y_continuous(name = "Renixx index", sec.axis = sec_axis(~.*5, name="Second Axis"))


#GEOPOLITICAL DATA
data_gpr_export <- read_excel("GitHub/PhD-Thesis/Poject paper 1/data_gpr_export.xls")
geo_i <- data_gpr_export %>% dplyr::select(1,5)
remove(data_gpr_export)
data <- geo_i[geo_i$month >= "2011-02-01" & geo_i$month <= "2023-01-01", ]
colnames(data)[1] <- "date"

data$loggeo <- log(data$GPRH)
data$test <- results$bsadf_panel
df_3 <- df_2[df_2$date >= "2011-02-01" & df_2$date <= "2023-01-01", ]
data$logstock <- df_3$logstock
ggplot(data, aes(date, loggeo))+geom_line()


#RATIO TEST
ratio<-data.frame(data$test,c95$`95%`)
ratio/ratio[,2]
data$ratio_test<-ratio$data.test


#CLIMATE CHANGE DATA
temp <- read.csv("~/GitHub/PhD-Thesis/Poject paper 1/Temperature anomalies.csv")
temp$Year<-seq(as.Date("1850-01-01"),as.Date("2023-01-01"),by='1 month')
data_2<-temp[temp$Year >= "2011-02-01" & temp$Year < "2023-01-01", ]
data$temp<-data_2$Value
ggplot(data, aes(date, temp))+geom_line()


# OIL price
#https://www.investing.com/commodities/crude-oil-historical-data
oil <- read.csv("~/GitHub/PhD-Thesis/Poject paper 1/Crude Oil WTI Futures Historical Data_daily.csv")
oil$Date <- mdy(oil$Date)
oil$logoil <- log(oil$Price)
ggplot(oil, aes(Date, Price))+geom_line()

oil_month <- oil %>% mutate(month = format(Date, "%m"), year = format(Date, "%Y")) %>% group_by(year,month) %>% summarize(Mean_m = mean(Price)) %>% as.data.frame
oil_m <- oil_month %>% mutate('Date' = make_date(year = year, month = month))
remove(oil_month)
ggplot(oil_m, aes(Date, Mean_m))+geom_line()

oil_m$log <- log(oil_m$Mean_m)
df_4 <- oil_m %>% dplyr::select(4,5)
crit_values_2 <- radf_mc_cv(nrow(df_4), nrep = 2000, seed=123)
results_2 <- radf(df_4)
datestamp(results_2, crit_values_2)
c95_2 <- as.data.frame(crit_values[["bsadf_cv"]])
autoplot(results_2, crit_values_2) + labs(title = "Oil GSADF")

autoplot(results_2,crit_values_2) + 
  labs(title = "GSADF and oil price", xlab= "Date", ylab="GSADF test") +
  geom_line(aes(x = Date, y = Mean_m/4),data = oil_m, color = 'black', inherit.aes = FALSE) +
  scale_y_continuous(name = "oil price", sec.axis = sec_axis(~.*2, name="Second Axis"))

oil_2<-oil_m[oil_m$Date >= "2011-02-01" & oil_m$Date < "2023-01-01", ]
data$logoil<-oil_2$log

#DIFF FOR VAR
data_df = cbind(diff(data$loggeo), diff(ratio$data.test), diff(data$temp))
colnames(data_df) <- c("geo_return", "test_return", 'temp_return')
plot.ts(data_df , main = "", xlab = "")
VARselect(data_df, type= "const", lag.max = 10)
fitvar1= VAR(data_df, p=1, type="both")
summary(fitvar1)

feir <- irf(fitvar1, impulse = "geo_return", response = c("test_return"), n.ahead = 8, ortho = FALSE, runs = 1000, seed = 003)
plot(feir)


#CATEGORICAL
renixx_m$dummy <- as.factor(ifelse(renixx_m$date >= "2011-08-01" & renixx_m$date <= "2013-01-01", 'NE',
                                   ifelse(renixx_m$date >= "2019-12-01" & renixx_m$date <= "2020-03-01"|              
                                            renixx_m$date >= "2020-07-01" & renixx_m$date <= "2021-05-01"|              
                                            renixx_m$date >= "2021-11-01" & renixx_m$date <= "2021-12-01", "PO", "ST")))


renixx_m2<-renixx_m[renixx_m$date>= "2011-02-01" & renixx_m$date <= "2023-01-01", ]
data$dummy=renixx_m2$dummy

data$geo_return<-(log(data$GPRH)-mean(log(data$GPRH)))
data$temp_return<-(log(data$temp)-mean(log(data$temp)))
data$oil_return<-(data$logoil-mean(data$logoil))
data$dummy_geo <- ifelse(data$geo_return >= 0, 1, 0)
data$dummy_temp <- ifelse(data$temp_return >= 0, 1, 0)
data$dummy_oil <- ifelse(data$oil_return >= 0, 1, 0)
data$dummy = relevel(data$dummy, ref = "ST")
model_1 <- multinom(dummy ~ dummy_geo + dummy_temp + dummy_oil, data = data)
summary(model_1)
tab_model(model_1)

table(data$dummy_oil,data$dummy)
table(data$dummy_geo,data$dummy)
table(data$dummy_temp,data$dummy)

#DUMMY
renixx_m$dummy <- as.factor(ifelse(renixx_m$date >= "2019-12-01" & renixx_m$date <= "2020-03-01"|              
                                     renixx_m$date >= "2020-07-01" & renixx_m$date <= "2021-05-01"|              
                                     renixx_m$date >= "2021-11-01" & renixx_m$date <= "2021-12-01", "1", "0"))

renixx_m2<-renixx_m[renixx_m$date>= "2011-02-01" & renixx_m$date <= "2023-01-01", ]
data$dummy=renixx_m2$dummy

model_2 <- glm(dummy ~ dummy_geo + dummy_temp + dummy_oil, data = data, family = binomial)
summary(model_2)
tab_model(model_2)

table(data$dummy, data$dummy_oil)

#OVERALL DATA
renixx$dummy <- as.factor(ifelse(renixx$date >= "2019-12-01" & renixx$date <= "2020-03-01"|              
                                     renixx$date >= "2020-07-01" & renixx$date <= "2021-05-01"|              
                                     renixx$date >= "2021-11-01" & renixx$date <= "2021-12-01", "1", "0"))

#renixx_2<-renixx[renixx$date>= "2011-02-01" & renixx$date <= "2023-01-01", ]

Week <- as.Date(cut(df$date, "week"))
renixx_w <- aggregate(close ~ Week, df, mean)

Week <- as.Date(cut(oil$Date, "week"))
oil_3 <- aggregate(Price ~ Week, oil, mean)

crit_values_w<-radf_mc_cv(nrow(renixx_w), seed = 0005)
results_w<-radf(renixx_w)
datestamp(results_w, crit_values_w)
c95<-as.data.frame(crit_values[["bsadf_cv"]])
autoplot(results_w, crit_values_w)
