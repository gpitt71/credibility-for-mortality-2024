rm(list=ls())
library(StMoMo)
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(demography)

# Full model ----

DK.mx <- hmd.mx(country = "ITA","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")

SWE.mx <- hmd.mx(country = "JPN","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")

# NOR.mx <- hmd.mx(country = "NOR","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")


ages.fit <- 50:89
years_fit_basic <- 1951:2003

series_code='male'

mortality_model_lc <- lc(link="log")
mortality_model_apc <- apc(link="log")
mortality_model_rh <- rh(link="log")

out <- NULL 

N_groups=2

forecasting_horizon <- 1
prediction_horizon <- 0:(15 - forecasting_horizon)
model_option <- "lc"
assign("mortality_model", get(paste0("mortality_model_", model_option)))
ix=0
years.fit <- min(years_fit_basic):max(years_fit_basic + ix)
E1 <- DK.mx$pop[[series_code]][as.character(ages.fit),as.character(years.fit)]
E2 <- SWE.mx$pop[[series_code]][as.character(ages.fit),as.character(years.fit)]
# E3 <- NOR.mx$pop[[series_code]][as.character(ages.fit),as.character(years.fit)]

D1 <- E1  * DK.mx$rate[[series_code]][as.character(ages.fit),as.character(years.fit)]
D2 <- E2  * SWE.mx$rate[[series_code]][as.character(ages.fit),as.character(years.fit)]
# D3 <- E3  * NOR.mx$rate[[series_code]][as.character(ages.fit),as.character(years.fit)]

Ext <- E1 #+ E2 + E3 
Dxt <-  D1 #+ D2 + D3 

list_of_extra_exposures <- list()

for (i in 2:N_groups) {
  
  tmp_list <- list()
  
  tmp_list[['Dxt']] <- get(paste0("D", i))
  tmp_list[['Ext']] <- get(paste0("E", i))
  list_of_extra_exposures[[i-1]] <- tmp_list
  
}


datahat <- structure(list(Dxt = Dxt, 
                          Ext = Ext, 
                          ages = ages.fit, 
                          years = years.fit, 
                          type = 'central', 
                          series = series_code, label = 'total'), 
                     class = "StMoMoData")



mortality_model_fit <- fit(mortality_model, 
                           data = datahat,
                           years.fit = years.fit,
                           ages.fit = ages.fit,
                           list_of_extra_exposures = list_of_extra_exposures)

year.predict <- max(years.fit)+forecasting_horizon

cv.arima.kt <- auto.arima(as.numeric(mortality_model_fit$kt), ic="bic")

if(model_option != "lc"){
  cv.arima.gc <- auto.arima(as.numeric(mortality_model_fit$gc), ic="bic")
  gc.order <- unname(arimaorder(cv.arima.gc))
}else{
  
  gc.order <- c(1,1,0)
}

mortality_model_forecast <- forecast(mortality_model_fit, 
                                     kt.method = "iarima",
                                     gc.order=gc.order,
                                     kt.order=unname(arimaorder(cv.arima.kt)),
                                     h = forecasting_horizon)


muxt_actual_1 <- DK.mx$rate[[series_code]][as.character(ages.fit),as.character(year.predict)] 
muxt_actual_2 <- SWE.mx$rate[[series_code]][as.character(ages.fit),as.character(year.predict)] 
# muxt_actual_3 <- NOR.mx$rate[[series_code]][as.character(ages.fit),as.character(year.predict)] 

if(forecasting_horizon>1){
  muxt_hat_predicted <- mortality_model_forecast$rates[as.character(ages.fit),as.character(year.predict)]
}else{
  muxt_hat_predicted <- mortality_model_forecast$rates
}



# Credibility model -----

ax_mx <- matrix(rep(mortality_model_fit$ax, dim(D1)[2]),
                                            nrow = dim(D1)[1],
                                            byrow= FALSE)

bx_mx <- matrix(rep(mortality_model_fit$bx, dim(D1)[2]),
                nrow = dim(D1)[1],
                byrow= FALSE)


kt_mx <- matrix(rep(mortality_model_fit$kt, dim(D1)[1]),
                nrow = dim(D1)[1],
                byrow= TRUE)

muxt_hat <- exp(ax_mx+bx_mx*kt_mx)

C1 <- apply(D1,1,sum)/apply(E1*muxt_hat,1,sum)

C2 <- apply(D2,1,sum)/apply(E2*muxt_hat,1,sum)

# C3 <- apply(D3,1,sum)/apply(E3*muxt_hat,1,sum)


C1_matrix <- matrix(rep(C1, dim(D1)[2]),
                    nrow = dim(D1)[1],
                    byrow= FALSE)

C2_matrix <- matrix(rep(C2, dim(D1)[2]),
                    nrow = dim(D1)[1],
                    byrow= FALSE)


# C3_matrix <- matrix(rep(C3, dim(D1)[2]),
#                     nrow = dim(D1)[1],
#                     byrow= FALSE)




lik_h0 <- sum(D1*log(E1*muxt_hat)-E1*muxt_hat-lfactorial(D1))+sum(D2*log(E2*muxt_hat)-E2*muxt_hat-lfactorial(D2))#+sum(D3*log(E3*muxt_hat)-E3*muxt_hat)
lik_h1 <- sum(D1*log(E1*muxt_hat*C1_matrix)-E1*muxt_hat*C1_matrix-lfactorial(D1))+sum(D2*log(E2*muxt_hat*C2_matrix)-E2*muxt_hat*C2_matrix-lfactorial(D2))#+sum(D3*log(E3*muxt_hat*C3_matrix)-E3*muxt_hat*C3_matrix)

lrt <- -2*(lik_h0- lik_h1)
lrt
qchisq(.999,df=N_groups*length(ages.fit))


tmp1 <- sum(D1*log(E1*muxt_hat)-E1*muxt_hat-lfactorial(D1))





