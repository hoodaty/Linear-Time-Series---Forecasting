datafile <- "C:/Users/DELL/OneDrive/Desktop/IPP Coursework M1/Time Series/data.csv"
data <- read.csv(datafile, header = FALSE,sep=",")

require(zoo)#practical and easy-to-use time series format (but larger)
require(tseries) #various functions on time series

#### Q1 ####

data[, 2] <- as.numeric(data[, 2])
df.source <- zoo(data[, 2])#convert the second column and name it as df.source into a time series
T <- length(df.source)#calculate the length of the time series
df <- df.source[1:T]#selecting the rows for plotting and storing them in df

####Plotting###
plot(df, xaxt="n", main = "Time Series Data") #represents df without the x-axis
# Define the sequence of dates from January 1990 to December 2018
dates <- seq(as.Date("1990-01-01"), by = "month", length.out = T)
# Set the x-axis ticks for year-wise plotting
axis(side = 1, at = seq(1, T, by = 12), labels = format(dates[seq(1, T, by = 12)], "%Y"))





###Q2####
d_df <- diff(df,1)#Differentiating to make it stationary
plot(cbind(df,d_df))
acf(d_df)


#decomposing to check trend and seasonality
df_s<- ts(df,frequency = 12) #converts a column from a data frame to a simple time series object.


# Decompose the time series
decomposed_ts <- decompose(df_s)
# Extract the trend component
trend_component <- decomposed_ts$trend
# Plot the trend component
plot(trend_component, main = "Trend Component of the Time Series", ylab = "Trend", xlab = "Time")








#checking stationarity using ADF Test

install.packages('fUnitRoots')
library('fUnitRoots')
require(fUnitRoots)





summary(lm(df ~ dates))
# ADF test
adf <- adfTest(df, lag=trunc((length(df)-1)^(1/3)), type="ct") 
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}


Qtests(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))


series <- df; kmax <- 348; adftype="ct"
adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests(adf@test$lm$residuals, 348, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}
adf <- adfTest_valid(df,348,adftype="ct")
Qtests(adf@test$lm$residuals, 48, fitdf = length(adf@test$lm$coefficients))
adf


# ADF test for differentiated series
summary(lm(d_df ~ dates[-1]))



adf1 <- adfTest(d_df, lag=trunc((length(d_df)-1)^(1/3)), type="nc")

Qtests(adf1@test$lm$residuals, 72, fitdf = length(adf1@test$lm$coefficients))

adf1 <- adfTest_valid(d_df,348,adftype="nc")
adf1



##pp test for stationarity
pp.test(d_df)
pp.test(df)

#KPSS Test
kpss.test(df,null='Trend')
kpss.test(df,null='Level')
kpss.test(d_df)

####Q3###
plot(cbind(df,d_df))








############PART-II####################



##########Q4########

#study of acf and pacf
x <- d_df
par(mfrow=c(1,2))
acf(d_df);pacf(d_df)




#choosing q=1 and p=3
pmax=3;qmax=1
mat <- matrix(NA,nrow=pmax+1,ncol=qmax+1) #empty matrix to fill
rownames(mat) <- paste0("p=",0:pmax) #renames lines
colnames(mat) <- paste0("q=",0:qmax) #renames columns
AICs <- mat #AIC matrix not filled non remplie
BICs <- mat #BIC matrix not filled non remplie
pqs <- expand.grid(0:pmax,0:qmax) #all possible combinations of p and q
for (row in 1:dim(pqs)[1]){ #loop for each (p,q)
  p <- pqs[row,1] #gets p
  q <- pqs[row,2] #gets q
  estim <- try(arima(x,c(p,0,q),include.mean = F)) #tries ARIMA estimation
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic #assigns the AIC
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) #assigns the BIC
}
AICs
AICs==min(AICs)
BICs
BICs==min(BICs)

#checking if well-adjusted
arima001 <- arima(d_df,c(0,0,1),include.mean=T)



#Checking Validity
Qtests(arima001$residuals,24,fitdf=3)

arima001

Box.test(arima001$residuals, lag=2, type="Ljung-Box", fitdf=1) 

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
Qtests(arima001$residuals, 24, 5) 


#Plotting residuals 
plot(arima001$residuals)

arima001 <- arima(d_df,c(0,0,1),include.mean=T)
arima001

###########PART-III#########

future_values <- forecast(arima001, h = 2,level=0.95)

plot(future_values)
future_values

# Print the future values
future_values
#confidence region

sigma2 <- arima001$sigma2
# Define the covariance matrix Σ
Sigma <- sigma2 * matrix(c(1, -0.5636, -0.5636, 1.31764496), nrow=2)

# Extract point forecasts and standard errors for X_{T+1} and X_{T+2}
point_forecasts <- future_values$mean[1:2]
se <- forecasted_values$se[1:2]

# Define the covariance matrix Σ
Sigma <- sigma2 * matrix(c(1, -0.5636, -0.5636, 1.31764496), nrow=2)

# Define the confidence level
alpha <- 0.05
chisq_val <- qchisq(1 - alpha, df=2)


# Compute the confidence ellipse centered at the forecasted values
ellipse_points <- ellipse(Sigma, centre=point_forecasts, level=1-alpha,t=chisq_val)

# Plotting the results
plot(ellipse_points, type='l', main="Time Series with Forecasts and Bivariate 95% Confidence Region",
     xlab="forecasted value of X_T+1|T",ylab="forecasted value of X_T+2|T",lwd=2)


abline(h=point_forecasts[2],v=point_forecasts[2],lwd=2)






##Test on Residuals


# Extract the residuals from the fitted ARIMA model
residuals <- residuals(arima001)


# Plot the residuals to visually inspect them
plot(residuals, main = "Residuals of ARIMA Model", ylab = "Residuals")


# Plot the histogram of the residuals
hist(residuals, breaks = 20, main = "Histogram of Residuals", xlab = "Residuals")


# Perform the Shapiro-Wilk test for normality
shapiro_test <- shapiro.test(residuals)
print(shapiro_test)

# Perform the Jarque-Bera test for normality
jarque_bera_test <- jarque.bera.test(residuals)
print(jarque_bera_test)

# Plot a Q-Q plot to visually check for normality
qqnorm(residuals)
qqline(residuals, col = "red")

# Check for causality of the solution
