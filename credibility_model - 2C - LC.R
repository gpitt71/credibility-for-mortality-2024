rm(list=ls())
library(StMoMo)
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(demography)

USA.mx <- hmd.mx(country = "USA","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")
USAdata <- StMoMoData(USA.mx,
                      series=names(USA.mx$rate)[2],
                      type = "central")


SWE.mx <- hmd.mx(country = "SWE","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")
SWEdata <- StMoMoData(SWE.mx,
                      series=names(SWE.mx$rate)[2],
                      type = "central")


ages.fit <- 55:89
years.fit <- 1951:2000

series='male'

E1 <- USA.mx$pop[[series]][as.character(ages.fit),as.character(years.fit)]
E2 <- SWE.mx$pop[[series]][as.character(ages.fit),as.character(years.fit)]
# E3 <- USA.mx$pop[[series]][as.character(ages.fit),as.character(years.fit)]
D1 <- E1  * USA.mx$rate[[series]][as.character(ages.fit),as.character(years.fit)]
D2 <- E2  * SWE.mx$rate[[series]][as.character(ages.fit),as.character(years.fit)]
# D3 <- E3  * USA.mx$rate[[series]][as.character(ages.fit),as.character(years.fit)]
Ext<-E1+E2#+E3
Dxt <-  D1 + D2 #+D3

datahat <- structure(list(Dxt = Dxt, 
                          Ext = Ext, 
                          ages = ages.fit, 
                          years = years.fit, 
                          type = 'central', 
                          series = series, label = 'total'), 
                     class = "StMoMoData")


LC <- lc(link="log")

LCModelfit <- fit(LC, 
             data = datahat,
             years.fit = years.fit,
             ages.fit = ages.fit)


#ages.fit is the same as ages.predict
years.predict <- 2001:2016
fcst.horizon <- length(years.predict)

LCModelfor <- forecast(LCModelfit, h = fcst.horizon)

E1.actual <- USA.mx$pop[[series]][as.character(ages.fit),as.character(years.predict)]
E2.actual <- SWE.mx$pop[[series]][as.character(ages.fit),as.character(years.predict)]
# E3.actual <- USA.mx$pop[[series]][as.character(ages.fit),as.character(years.predict)]

Ext.actual <-E1.actual+E2.actual#+E3.actual

Dxt.actual <- E1.actual * USA.mx$rate[[series]][as.character(ages.fit),as.character(years.predict)] + E2.actual  * SWE.mx$rate[[series]][as.character(ages.fit),as.character(years.predict)]#+ E3.actual  * USA.mx$rate[[series]][as.character(ages.fit),as.character(years.predict)]
Dxt.hat <-  Ext.actual*LCModelfor$rates

mse0 <- mean((Dxt.hat-Dxt.actual)^2)


# Country one model

weights_1_dot <- apply(E1,MARGIN = 1,sum)

F_1 <- apply(D1*E1/E1,MARGIN = 1,sum)/weights_1_dot

# This is equivalent to: 
# LCModelfit1 <- fit(LC, 
#               data = USAdata,
#               years.fit = years.fit,
#               ages.fit = ages.fit)
# F_1==exp(LCModelfit1$ax)


F_1_mx <- matrix(rep(F_1,fcst.horizon),byrow = F,nrow = length(ages.fit))

# Country two model

weights_2_dot <- apply(E2,MARGIN = 1,sum)

F_2 <- apply(E2*D2/E2,MARGIN = 1,sum)/weights_2_dot

F_2_mx <- matrix(rep(F_2,fcst.horizon),byrow = F,nrow = length(ages.fit))

weights_dot_dot <- weights_1_dot + weights_2_dot

# estimation of structural parameters 

I = 2

lambda0 <- LCModelfor$fitted[,1]

#equivalent to: F_1*tmp1+F_2*tmp2+F_3*tmp3

tmp1 <- weights_1_dot/weights_dot_dot
tmp2 <- weights_2_dot/weights_dot_dot
c.den <- tmp1*(1-tmp1)+tmp2*(1-tmp2)

c.init <- ((I-1)/I)*(1/c.den)

T.init <- (I/(I-1))*(tmp1*(F_1-lambda0)^2+tmp2*(F_2-lambda0)^2)

tau0 <- c.init*(T.init-(I*lambda0/weights_dot_dot))

max.iter=100
l1 <- list()
l2 <- list()
l3 <- list()
for(iter in 1:max.iter){
  
  k.iter <- lambda0/tau0
  alpha1.iter <- weights_1_dot/(weights_1_dot+k.iter)
  alpha2.iter <- weights_2_dot/(weights_2_dot+k.iter)

  l1[[iter]] <- k.iter
  l2[[iter]] <- alpha1.iter
  l3[[iter]] <- alpha2.iter
  
  alpha.dot <- alpha1.iter+alpha2.iter
  # lambda0 <- (alpha1.iter/alpha.dot)*F_1 + (alpha2.iter/alpha.dot)*F_2
  tau0 <- max(c.init*(T.init-((I*lambda0)/weights_dot_dot)),0)
  
}

alpha1.iter <-  weights_1_dot/(weights_1_dot+k.iter)
alpha2.iter <-  weights_2_dot/(weights_2_dot+k.iter)

alpha1.mx <- matrix(rep(alpha1.iter,fcst.horizon),byrow = F,nrow = length(ages.fit))
alpha2.mx <- matrix(rep(alpha2.iter,fcst.horizon),byrow = F,nrow = length(ages.fit))

credible.mu.hat.1 <- LCModelfor$rates+alpha1.mx*(F_1_mx-LCModelfor$rates)
credible.mu.hat.2 <- LCModelfor$rates+alpha2.mx*(F_2_mx-LCModelfor$rates)

credible.Dxt.hat <- E1.actual*credible.mu.hat.1 + E2.actual*credible.mu.hat.2

mse1 <- mean((credible.Dxt.hat-Dxt.actual)^2)

##

LCModelfit1 <- fit(LC,
              data = USAdata,
              years.fit = years.fit,
              ages.fit = ages.fit)

LCModelfit2 <- fit(LC,
              data = SWEdata,
              years.fit = years.fit,
              ages.fit = ages.fit)

LCModelfor1 <- forecast(LCModelfit1, h = fcst.horizon)
LCModelfor2 <- forecast(LCModelfit2, h = fcst.horizon)


Dxt.hat.sep <-  E1.actual*LCModelfor1$rates+E2.actual*LCModelfor2$rates
mse2 <- mean((Dxt.hat.sep-Dxt.actual)^2)

print(c(mse0,mse1,mse2))

