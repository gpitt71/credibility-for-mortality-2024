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


lc.model = lc(link="log")

LCModelfit <- fit(lc.model, 
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


## -----

LCModelfit1 <- fit(lc.model,
                   data = USAdata,
                   years.fit = years.fit,
                   ages.fit = ages.fit)

LCModelfit2 <- fit(lc.model,
                   data = SWEdata,
                   years.fit = years.fit,
                   ages.fit = ages.fit)

LCModelfor1 <- forecast(LCModelfit1, h = fcst.horizon)
LCModelfor2 <- forecast(LCModelfit2, h = fcst.horizon)

I = 2

T.init <- (I/(I-1))*((E1.actual/(E1.actual+E2.actual))*(LCModelfor1$rates-LCModelfor$rates)^2+(E2.actual/(E1.actual+E2.actual))*(LCModelfor2$rates-LCModelfor$rates)^2)/2

tmp1=E1.actual/(E1.actual+E2.actual)
tmp2=E2.actual/(E1.actual+E2.actual)
                
c.den <- tmp1*(1-tmp1)+tmp2*(1-tmp2)

c.init <- ((I-1)/I)*(1/c.den)

tau0 <- c.init*(T.init-(I*LCModelfor$rates/(E1.actual+E2.actual)))

k=apply(LCModelfor$rates/(tau0),MARGIN=2,pmax,0)

alpha1.mx <- E1.actual/(E1.actual+k)
alpha2.mx <- E2.actual/(E2.actual+k)

credible.mu.hat.1 <- LCModelfor$rates+alpha1.mx*(LCModelfor1$rates-LCModelfor$rates)
credible.mu.hat.2 <- LCModelfor$rates+alpha2.mx*(LCModelfor2$rates-LCModelfor$rates)

credible.Dxt.hat <- E1.actual*credible.mu.hat.1 + E2.actual*credible.mu.hat.2

mse1 <- mean((credible.Dxt.hat-Dxt.actual)^2)


Dxt.hat.sep <-  E1.actual*LCModelfor1$rates+E2.actual*LCModelfor2$rates
mse2 <- mean((Dxt.hat.sep-Dxt.actual)^2)

print(c(mse0,mse1,mse2))








