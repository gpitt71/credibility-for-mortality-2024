rm(list=ls())
library(StMoMo)
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(demography)

USA.mx <- hmd.mx(country = "USA","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")
USAdata <- StMoMoData(USA.mx,
                      series=names(USA.mx$rate)[1],
                      type = "central")


SWE.mx <- hmd.mx(country = "SWE","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")
SWEdata <- StMoMoData(SWE.mx,
                      series=names(SWE.mx$rate)[1],
                      type = "central")

ITA.mx <- hmd.mx(country = "ITA","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")
ITAdata <- StMoMoData(ITA.mx,
                     series=names(USA.mx$rate)[1],
                     type = "central")

ages.fit <- 55:89
years.fit <- 1951:2000

series='total'

E1 <- USA.mx$pop[[series]][as.character(ages.fit),as.character(years.fit)]
E2 <- SWE.mx$pop[[series]][as.character(ages.fit),as.character(years.fit)]
E3 <- ITA.mx$pop[[series]][as.character(ages.fit),as.character(years.fit)]
D1 <- E1  * USA.mx$rate[[series]][as.character(ages.fit),as.character(years.fit)]
D2 <- E2  * SWE.mx$rate[[series]][as.character(ages.fit),as.character(years.fit)]
D3 <- E3  * ITA.mx$rate[[series]][as.character(ages.fit),as.character(years.fit)]
Ext<-E1+E2+E3
Dxt <-  D1 + D2 +D3

datahat <- structure(list(Dxt = Dxt, 
                          Ext = Ext, 
                          ages = ages.fit, 
                          years = years.fit, 
                          type = 'central', 
                          series = series, label = 'total'), 
                     class = "StMoMoData")


LC <- lc(link = "log")
LCfit <- fit(LC, 
             data = datahat,
             years.fit = years.fit,
             ages.fit = ages.fit)



#ages.fit is the same as ages.predict
years.predict <- 2001:2016
fcst.horizon <- length(years.predict)

LCfor <- forecast(LCfit, h = fcst.horizon)

E1.actual <- USA.mx$pop[[series]][as.character(ages.fit),as.character(years.predict)]
E2.actual <- SWE.mx$pop[[series]][as.character(ages.fit),as.character(years.predict)]
E3.actual <- ITA.mx$pop[[series]][as.character(ages.fit),as.character(years.predict)]

Ext.actual <-E1.actual+E2.actual+E3.actual

Dxt.actual <- E1.actual * USA.mx$rate[[series]][as.character(ages.fit),as.character(years.predict)] + E2.actual  * SWE.mx$rate[[series]][as.character(ages.fit),as.character(years.predict)]+ E3.actual  * ITA.mx$rate[[series]][as.character(ages.fit),as.character(years.predict)]
Dxt.hat <-  Ext.actual*LCfor$rates

mse0 <- mean((Dxt.hat-Dxt.actual)^2)


# Country one model

weights_1_dot <- apply(E1,MARGIN = 1,sum)

F_1 <- apply(D1*E1/E1,MARGIN = 1,sum)/weights_1_dot

F_1_mx <- matrix(rep(F_1,fcst.horizon),byrow = F,nrow = length(ages.fit))

# Country two model

weights_2_dot <- apply(E2,MARGIN = 1,sum)

F_2 <- apply(D2*E2/E2,MARGIN = 1,sum)/weights_2_dot

F_2_mx <- matrix(rep(F_2,fcst.horizon),byrow = F,nrow = length(ages.fit))

# Country three model

weights_3_dot <- apply(E3,MARGIN = 1,sum)

weights_dot_dot <- weights_1_dot + weights_2_dot + weights_3_dot

F_3 <- apply(D3*E3/E3,MARGIN = 1,sum)/weights_3_dot

F_3_mx <- matrix(rep(F_3,fcst.horizon),byrow = F,nrow = length(ages.fit))

# estimation of structural parameters 

I = 3

tmp1 <- weights_1_dot/weights_dot_dot
tmp2 <- weights_2_dot/weights_dot_dot
tmp3 <- weights_3_dot/weights_dot_dot
c.den <- tmp1*(1-tmp1)+tmp2*(1-tmp2)+tmp3*(1-tmp3)

c.init <- ((I-1)/I)*(1/c.den)

F_bar <- LCfor$fitted #F_1*tmp1+F_2*tmp2+F_3*tmp3

T.init <- (I/(I-1))*(tmp1*(F_1-F_bar)^2+tmp2*(F_2-F_bar)^2+tmp3*(F_3-F_bar)^2)

max.iter = 100

lambda.0 <- F_bar
tau2.0 <- c.init*(T.init-((I*lambda.0)/weights_dot_dot))

l1 <- list()
l2 <- list()
l3 <- list()
l4 <- list()

for(iter in 1:max.iter){
  
  k.iter <- lambda.0/tau2.0
  alpha1.iter <- weights_1_dot/(weights_1_dot+k.iter)
  alpha2.iter <- weights_2_dot/(weights_2_dot+k.iter)
  alpha3.iter <- weights_3_dot/(weights_3_dot+k.iter)
  
  l1[[iter]] <- k.iter
  l2[[iter]] <- alpha1.iter
  l3[[iter]] <- alpha2.iter
  l4[[iter]] <- alpha3.iter
  
  
  alpha.dot <- alpha1.iter+alpha2.iter+alpha3.iter
  lambda.0 <- (alpha1.iter/alpha.dot)*F_1 + (alpha2.iter/alpha.dot)*F_2+ (alpha3.iter/alpha.dot)*F_3
  tau2.0 <- c.init*(T.init-((I*lambda.0)/weights_dot_dot))
}

alpha1.mx <- matrix(rep(alpha1.iter,fcst.horizon),byrow = F,nrow = length(ages.fit))
alpha2.mx <- matrix(rep(alpha2.iter,fcst.horizon),byrow = F,nrow = length(ages.fit))
alpha3.mx <- matrix(rep(alpha3.iter,fcst.horizon),byrow = F,nrow = length(ages.fit))

credible.mu.hat.1 <- LCfor$rates+alpha1.mx*(F_1_mx-LCfor$rates)
credible.mu.hat.2 <- LCfor$rates+alpha2.mx*(F_2_mx-LCfor$rates)
credible.mu.hat.3 <- LCfor$rates+alpha3.mx*(F_3_mx-LCfor$rates)

credible.Dxt.hat <- E1.actual*credible.mu.hat.1 + E2.actual*credible.mu.hat.2 + E3.actual*credible.mu.hat.3

mse1 <- mean((credible.Dxt.hat-Dxt.actual)^2)


print(c(mse0,mse1)/1e05)




