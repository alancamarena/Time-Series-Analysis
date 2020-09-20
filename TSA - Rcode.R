setwd("C:/Users/Alan Camarena/Desktop")  # set working directory
maryland = read.table("drug induced deaths~maryland.csv", sep = ",", header = FALSE, skip = 1, nrows = 204)   # read data 
mary <- ts(maryland[,2], start = c(1999,1), frequency = 12) #convert the data into a time series 
ts.plot(mary, ylab = " Deaths per month", xlab = "Year date",xlim=c(1999,2016),  main = "Drug induced deaths per month (State of Maryland)") #time series plot

op <- par(mfrow = c(1,2)) #keep it together
acf(mary, lag.max = 144, main = " ")
pacf(mary, lag.max = 144, main = " ")
title ("Deaths sample ACF/PACF", outer = TRUE, line = -2 ) # center title
par(op) # tie the knot


require(MASS) #BoxCox test  
bcTransform <- boxcox(mary ~ as.numeric(1:length(mary)))
lambda = bcTransform$x[which(bcTransform$y == max(bcTransform$y))]
# -.7878788 (our max lambda), I consider using the optimal lambda or the recip

y.max <- (1/lambda)*(mary^lambda - 1) #box cox
ts.plot(y.max, main = "box-cox") #variance of 6.61x-5

y.rec = mary^(-1)  
ts.plot(y.rec, main = "Reciprocal of -1") #variance of 1.11x-5

y.ori = mary^(lambda)
ts.plot(y.ori, main = "power transform of max lambda") # variance is 4.1x-5

y1 = diff(y.max, 1) #difference at lag  = 1 to remove trend component 
ts.plot(y1, main = "De- trended and Transformed T.S.", ylab = " ")
abline(h = mean(y1), lty = 2, col = "red")
#variance of y.max is .0000661
#var(y1) # variance goes down 
#mean(y1) mean almost 0
#difference one more time to see if variance goes up 
#var(y1) decreases
y2 = diff(y1, 1) # differencing again to check variance 
#var(y2) # variance increases, 1 difference is sufficient 

#might want to leave big enough to observe, adjust size
acf(y1, lag.max = 36, main = "Differenced & Transformed ACF")
pacf(y1, lag.max = 36, main = "Differenced & Transformed PACF")

#ar(y1, aic = TRUE,  method = "mle") # is invertible and casual
#spec.arma(ar = c(-.7154,-.5872,-.4882,-.3438,-.4637,-.3811,-.3286,-.2626,-.3495,-.2770,-.1559 )) # no errors returned 
#ar(y1,aic = TRUE, method = "yule-walker") # is invertible and casual 
#spec.arma(ar = c(-.7108,-.5810,-.4807,-.3278,-.4484,-.3637,-.3151,-.2498,-.3364,-.2646,-.1494 )) # no errors returned 

library(qpcR)
fit_ar11 = arima(y1, order = c(11,0,0), method = "ML")
fit_ma9 = arima(y1, order = c(0,0,9), method = "ML")
#AICc(fit_ar11) 
#AICc(fit_ma9)   

aiccs <- matrix(NA, nr = 10, nc = 10) 
dimnames(aiccs) = list(p=0:9, q=0:9) 
for(p in 0:9) 
{
  for(q in 0:9) 
  {
    aiccs[p+1,q+1] = AICc(arima(y1, order = c(p,0,q), method="ML")) 
  }
}
aiccs

#r work to fit models 
fit35 = arima(y.max, order= c(3,1,5), method = "ML", xreg = 1:length(y.max) )#casual and invertible 
fit13 = arima(y.max, order = c(1,1,3), method = "ML")#casual and invertible 
fit11 = arima(y.max, order = c(1,1,1), method = "ML")

#code for zeroing out coefficients for each model 
#ARIMA(3,5)
#check every confidence interval individually 
c(-.8136 - 1.96*.0930, -.8136 + 1.96*.0930) #1  good      (AR PART)
c(-.8489 - 1.96*.0555, -.8489 + 1.96*.0555) #2  good
c(-.8208 - 1.96*.0876, -.8208 + 1.96*.0876) #3  good, all ar part is ok
c( .1066 - 1.96*.1092,  .1066 + 1.96*.1092) #1, zero out   (MA PART)
c( .2348 - 1.96*.0984,  .2348 + 1.96*.0984) #2, good
c( .1382 - 1.96*.1091,  .1382 + 1.96*.1091) #3, zero out
c(-.6642 - 1.96*.0854, -.6642 + 1.96*.0854) #4  good
c(-.2763 - 1.96*.0722, -.2763 + 1.96*.0722) #5  good

#ARIMA(1,3)
#check confidence interval to zero out coefficients
c( -.9132- 1.96*.0716, -.9132 + 1.96*.0716) 
c( .2194- 1.96*.0984, .2194 + 1.96*.0984) 
c( -.7379- 1.96*.0738, -.7379 + 1.96*.0738)
c( -.1956 - 1.96*.0726, -.1953 + 1.96*.0726) #all coefficients are ok to use


#ARIMA(3,5) 
ts.plot(residuals(fit35), main = "fitted residuals for ARIMA(3,1,5)") # the residuals look like white noise, stationary 
abline(h = mean(residuals(fit35)), lty = 2, col = "red")
acf(residuals(fit35), lag.max = 36, main = "ACF~ARIMA(3,1,5)")
pacf(residuals(fit35), lag.max = 36, main = "PACF~ARIMA(3,1,5)")

#construct new arma/arima model 
j33 = arima(y.max, order = c(3,1,5), seasonal = list(order = c(0,0,0), period = 12), fixed = c(NA,NA,NA,0,NA,0,NA,NA)) #need to zero out some coefficients pheta 1, pheta 3
#j33 # model's AICc does improve

#attempt diagnostic checking again. Residuals plot, acf/pacf
ts.plot(residuals(j33), main = "fitted residuals ~ARIMA(3,1,3)") # the residuals look like white noise, stationary 
abline(h = mean(residuals(j33)), lty = 2, col = "red")
acf(residuals(j33), lag.max = 48, main = "ACF~ ARIMA(3,1,3)")
pacf(residuals(j33), lag.max = 48, main = " PACF ~ ARIMA(3,1,3)") # there is some improvement in that coefficient 

hist((residuals(fit13)), prob = TRUE, main = "Hisogram Normality check") # distribution looks normal
lines(density(residuals(fit13), adjust = 2), col = "Red", lwd = 2) # histogram follows normality 

#spec.arma(ar = c(-0.6933, -0.8712, -0.7052), ma = c(0.3720 , -0.5467, -0.2792)) #it is invertible and casual(no warning message)
plot.roots(NULL, polyroot(c(1 , 0.6933 , 0.8712 , 0.7052 )), main = "Roots of AR part")
plot.roots(NULL, polyroot(c(1,0.3720, -0.5467, -0.2792)), main = "Roots of MA part")
#passes every test except shapiro, my normalitty 
Box.test(residuals(j33), lag = 14, type =c("Box-Pierce"), fitdf = 6) 
Box.test(residuals(j33), lag = 14, type = c("Ljung-Box"), fitdf = 6)
Box.test(residuals(j33)^2, lag = 14, type =c("Ljung-Box"), fitdf = 0) #check for linearity
#doesn't pass
shapiro.test(residuals(j33))

#ARMA(1,3) 

ts.plot(residuals(fit13), main = "fitted residuals ~ ARIMA(1,1,3)")  
abline(h = mean(residuals(fit13)), lty = 2, col = "red")
acf(residuals(fit13), lag.max = 48, main = "ACF~ARIMA(1,1,3)")
pacf(residuals(fit13), lag.max = 48, main = "PACF~ARIMA(1,1,3) ")

hist((residuals(fit13)), prob = TRUE, main = "Hisogram Normality check") # distribution looks normal
lines(density(residuals(fit13), adjust = 2), col = "Red", lwd = 2) # histogram follows normality 


#spec.arma(ar = c(-0.9117), ma = c(0.2148 , -0.7408, -0.1982)) # returns error message if not invertible and casual 
plot.roots(NULL, polyroot(c(1,0.9132)), main = "Roots of AR part")
plot.roots(NULL, polyroot(c(1,0.2194 , -0.7379, -0.1953)), main = "Roots of MA part") #plot.roots used for invertibility and casual checking

#passes every test except shapiro,  normality 
Box.test(residuals(fit13), lag = 14, type =c("Box-Pierce"), fitdf = 4) 
Box.test(residuals(fit13), lag = 14, type = c("Ljung-Box"), fitdf = 4)
Box.test(residuals(fit13)^2, lag = 14, type =c("Ljung-Box"), fitdf = 0) #check for linearity
#doesn't pass
shapiro.test(residuals(fit13))

#arima(1,1,1) 

ts.plot(residuals(fit11), main = "fitted residuals~ ARIMA(1,1,1)")  #looks like WN
abline(h = mean(residuals(fit11)), lty = 2, col = "red")
acf(residuals(fit11), lag.max = 48, main = "ACF~ ARIMA(1,1,1) ") # 3 residuals are out of CI
pacf(residuals(fit11), lag.max = 48, main = "PACF~ ARIMA(1,1,1) ")# 3 residuals are out of CI
hist((residuals(fit11)), prob = TRUE, main = "Hisogram Normality check") # distribution looks normal
lines(density(residuals(fit11), adjust = 2), col = "Red", lwd = 2) # histogram follows normality 


#spec.arma(ar = c(0.1359), ma = c(-0.8474)) # returns error message if not invertible and casual 
plot.roots(NULL, polyroot(c(1,-0.1359)), main = "Roots of AR part") #is invertible and casual 
plot.roots(NULL, polyroot(c(1,-0.8474)), main = "Roots of MA part") 

#again, passes everything except for shapiro wilk 
Box.test(residuals(fit11), lag = 14, type =c("Box-Pierce"), fitdf = 2) 
Box.test(residuals(fit11), lag = 14, type = c("Ljung-Box"), fitdf = 2)
Box.test(residuals(fit11)^2, lag = 14, type =c("Ljung-Box"), fitdf = 0) 
shapiro.test(residuals(fit11))
#read table for to see if our forecasting model is accurat, thi section entails  
maryland2 = read.table("take away 12 observations.csv", sep = ",", header = FALSE, skip = 1, nrows = 192)   # read data 
maryshort <- ts(maryland2[,2], start = c(1999,1), frequency = 12) #convert the data into a time series 
ts.plot(maryshort)
bcTransform1 <- boxcox(maryshort ~ as.numeric(1:length(maryshort)))
lambda1 = bcTransform$x[which(bcTransform$y == max(bcTransform$y))]
#transformation 
maryshortlog = log(maryshort)
ts.plot(maryshort) 
#difference
ymary = diff(maryshortlog, 1) #difference at lag  = 1 to remove trend component 
ts.plot(ymary, main = "De- trended and Transformed T.S.", ylab = " ")
abline(h = mean(y1), lty = 2, col = "red")
short35 = arima(maryshortlog, order= c(3,1,5), method = "ML", xreg = 1:length(maryshortlog) )
short35
aiccs1 <- matrix(NA, nr = 10, nc = 10) 
dimnames(aiccs1) = list(p=0:9, q=0:9) 
for(p in 0:9) 
{
  for(q in 0:9) 
  {
    aiccs1[p+1,q+1] = AICc(arima(ymary, order = c(p,0,q), method="ML")) 
  }
}
aiccs1
#ARIMA(3,1,5) is also adequate for model with 12 less observations
#FORECASTING, lets try to forecast 1 year
pred.tr = predict(fit35, n.ahead = 12,newxreg =  (length(y.max)+1):(length(y.max)+12)) #newxreg=(length(y.log)+1):length(y.log) + 10)
u.tr = pred.tr$pred + 2*pred.tr$se #upper confidence interval
l.tr = pred.tr$pred - 2*pred.tr$se #lower confidence interval 
ts.plot(y.max, xlim = c(1999,2019), ylim = c(min(y.max),max(u.tr)), main = "Transformed data forecasting")
lines(u.tr, col = "blue", lty = "dashed")# CI prediciton lines
lines(l.tr, col = "blue", lty = "dashed")
lines(pred.tr$pred, col = "red", lty = "dashed") #main prediciton line 

#attempt to produce original data
predict.original =  (lambda*(pred.tr$pred) + 1)^(1/lambda) 
U =  (lambda*(u.tr) + 1)^(1/lambda)
L =  (lambda*(l.tr) + 1)^(1/lambda)
ts.plot(mary,xlim = c(1999,2017), ylim = c(min(mary),max(U)), main = "Original data forecasting")
lines(U, col = "blue", lty = 2)
lines(L, col = "blue", lty = 2)
lines(predict.original, col = "red", pch = 16, lty = "dashed")
#now I will attempt to forecast 2015 of my original data ~ short information 
pred.mr = predict(short35, n.ahead = 12,newxreg =  (length(maryshortlog)+1):(length(maryshortlog)+12))
#newxreg=(length(y.log)+1):length(y.log) + 10)
u.tr1 = pred.mr$pred + 2*pred.mr$se
l.tr1 = pred.mr$pred - 2*pred.mr$se
ts.plot(maryshortlog, xlim = c(1999,2019), ylim = c(min(maryshortlog),max(u.tr1)))
lines(u.tr1, col = "blue", lty = "dashed")
lines(l.tr1, col = "blue", lty = "dashed")
lines(pred.mr$pred, col = "red", lty = "dashed")

#attempt to produce original data
predict.original1 =  exp(pred.mr$pred)
U1 =  exp(u.tr1)
L1 =  exp(l.tr1)
ts.plot(maryshort,xlim = c(1999,2017), ylim = c(min(maryshort),max(U1)))
lines(U1, col = "blue", lty = 2)
lines(L1, col = "blue", lty = 2)
lines(predict.original1, col = "red", pch = 16, lty = "dashed")
#original time series plot for comparison 
ts.plot(mary, ylab = " Deaths per month", xlab = "Year date",xlim=c(1999,2016),  main = "Drug induced deaths per month (State of Maryland)")