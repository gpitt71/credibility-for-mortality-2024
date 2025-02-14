rm(list=ls())
library(StMoMo)
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(demography)

# Full model ----
series_code='female'

DK.mx <- hmd.mx(country = "DNK","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")
D1data <- StMoMoData(DK.mx,
                     series=series_code,
                     type = "central")

SWE.mx <- hmd.mx(country = "SWE","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")
D2data <- StMoMoData(SWE.mx,
                     series=series_code,
                     type = "central")

NOR.mx <- hmd.mx(country = "NOR","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")
D3data <- StMoMoData(NOR.mx,
                     series=series_code,
                     type = "central")

FIN.mx <- hmd.mx(country = "FIN","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")
D4data <- StMoMoData(FIN.mx,
                     series=series_code,
                     type = "central")

ICE.mx <- hmd.mx(country = "ISL","gabriele.pittarello@uniroma1.it", "pm869fqW8nxozS8!")
D5data <- StMoMoData(ICE.mx,
                     series=series_code,
                     type = "central")

ages.fit <- 50:89
years_fit_basic <- 1951:2003



mortality_model_lc <- lc(link="log")
mortality_model_apc <- apc(link="log")
mortality_model_rh <- rh(link="log")

out <- NULL 

N_groups=3

for(forecasting_horizon in c(1,5,12)){
  
  prediction_horizon <- 0:(15-forecasting_horizon)
  
  for(model_option in c("lc","apc","rh")){
    
    mse_0 <- mse_1 <- mse_2 <- NULL
    
    assign("mortality_model",get(paste0("mortality_model_",model_option)))
    
    for(ix in prediction_horizon){
      
      years.fit<- min(years_fit_basic):max(years_fit_basic+ix)
      
      E1 <- DK.mx$pop[[series_code]][as.character(ages.fit),as.character(years.fit)]
      E2 <- SWE.mx$pop[[series_code]][as.character(ages.fit),as.character(years.fit)]
      E3 <- NOR.mx$pop[[series_code]][as.character(ages.fit),as.character(years.fit)]
      E4 <- FIN.mx$pop[[series_code]][as.character(ages.fit),as.character(years.fit)]
      E5 <- ICE.mx$pop[[series_code]][as.character(ages.fit),as.character(years.fit)]
      
      
      
      D1 <- E1  * DK.mx$rate[[series_code]][as.character(ages.fit),as.character(years.fit)]
      D2 <- E2  * SWE.mx$rate[[series_code]][as.character(ages.fit),as.character(years.fit)]
      D3 <- E3  * NOR.mx$rate[[series_code]][as.character(ages.fit),as.character(years.fit)]
      D4 <- E4  * FIN.mx$rate[[series_code]][as.character(ages.fit),as.character(years.fit)]
      D5 <- E5  * ICE.mx$rate[[series_code]][as.character(ages.fit),as.character(years.fit)]
      
            
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
      muxt_actual_3 <- NOR.mx$rate[[series_code]][as.character(ages.fit),as.character(year.predict)] 
      
      if(forecasting_horizon>1){
        muxt_hat_predicted <- mortality_model_forecast$rates[as.character(ages.fit),as.character(year.predict)]
      }else{
        muxt_hat_predicted <- mortality_model_forecast$rates
      }
      
      
      mse_0_ix <- 0
      
      for (i in 1:N_groups) {
        mse_0_ix <- mse_0_ix + (abs((get(paste0("muxt_actual_", i))-muxt_hat_predicted)))#/get(paste0("muxt_actual_", i))
        
      }
      
      mse_0 <- c(mse_0,sum(mse_0_ix))
      
      # Credibility model -----
      
      muxt_hat <- mortality_model_forecast$fitted
      
      C1 <- apply(D1,1,sum)/apply(E1*muxt_hat,1,sum)
      
      C2 <- apply(D2,1,sum)/apply(E2*muxt_hat,1,sum)
      
      C3 <- apply(D3,1,sum)/apply(E3*muxt_hat,1,sum)
      
      C4 <- apply(D4,1,sum)/apply(E4*muxt_hat,1,sum)
      
      C5 <- apply(D5,1,sum)/apply(E5*muxt_hat,1,sum)
      
      Fxt_1 <- D1/E1
      Fxt_2 <- D2/E2
      Fxt_3 <- D3/E3
      Fxt_4 <- D4/E4
      Fxt_5 <- D5/E5
      
      varthetax_1 <- apply(((Fxt_1-muxt_hat)^2)/(muxt_hat^2),1,mean)
      varthetax_2 <- apply(((Fxt_2-muxt_hat)^2)/(muxt_hat^2),1,mean)
      varthetax_3 <- apply(((Fxt_3-muxt_hat)^2)/(muxt_hat^2),1,mean)
      varthetax_4 <- apply(((Fxt_4-muxt_hat)^2)/(muxt_hat^2),1,mean)
      varthetax_5 <- apply(((Fxt_5-muxt_hat)^2)/(muxt_hat^2),1,mean)
      
      Z_1 <- 1/(1+varthetax_1*apply(E1*muxt_hat,1,sum))
      
      Z_2 <- 1/(1+varthetax_2*apply(E2*muxt_hat,1,sum))
      
      Z_3 <- 1/(1+varthetax_3*apply(E3*muxt_hat,1,sum))
      
      Z_4 <- 1/(1+varthetax_4*apply(E4*muxt_hat,1,sum))
      
      Z_5 <- 1/(1+varthetax_5*apply(E5*muxt_hat,1,sum))
      
      muhat1 <- Z_1*muxt_hat_predicted + (1-Z_1)* C1 * muxt_hat_predicted
      muhat2 <- Z_2*muxt_hat_predicted + (1-Z_2)* C2 * muxt_hat_predicted
      muhat3 <- Z_3*muxt_hat_predicted + (1-Z_3)* C3 * muxt_hat_predicted
      muhat4 <- Z_4*muxt_hat_predicted + (1-Z_4)* C4 * muxt_hat_predicted
      muhat5 <- Z_5*muxt_hat_predicted + (1-Z_5)* C5 * muxt_hat_predicted
      
      mse_1_ix <- 0
      
      for (i in 1:N_groups) {
        
        mse_1_ix <- mse_1_ix + (abs(get(paste0("muxt_actual_", i))-get(paste0("muhat", i))))#/get(paste0("muxt_actual_", i))
        
      }
      
      mse_1 <- c(mse_1,sum(mse_1_ix))
      
      # Separated model ----
      
      mse_2_ix <- 0
      
      for(i in 1:N_groups){
        
        popspecdata <- get(paste0("D",i,"data"))
        
        assign(paste0("mortality_model_fit",i),
               fit(mortality_model,
                   data = popspecdata,
                   years.fit = years.fit,
                   ages.fit = ages.fit))
        
        assign(paste0("cv.arima.kt_", i), auto.arima(as.numeric(get(paste0("mortality_model_fit",i))$kt)))
        
        
        if(model_option != "lc"){
          assign(paste0("cv.arima.gc_", i), auto.arima(as.numeric(get(paste0("mortality_model_fit",i))$gc)))
          assign(paste0("gc.order_", i), unname(arimaorder(get(paste0("cv.arima.gc_", i)))))
          
        }else{
          assign(paste0("gc.order_", i),c(1,1,0))
        }
        
        
        mortality_model_fit
        
        error_check <- tryCatch(assign(paste0("mortality_model_forecast",i) , forecast(get(paste0("mortality_model_fit",i)), 
                                                                                       gc.order = get(paste0("gc.order_",i)),
                                                                                       kt.method = "iarima",
                                                                                       kt.order=unname(arimaorder(get(paste0("cv.arima.kt_",i)))),
                                                                                       h = forecasting_horizon)),
                                error = function(e){
                                  
                                  return("mrw")
                                  
                                  
                                })
        
        
        if(is.character(error_check)){
          
          assign(paste0("mortality_model_forecast",i) , 
                 forecast(get(paste0("mortality_model_fit",i)),
                          h = forecasting_horizon))
        }else{
          
          assign(paste0("mortality_model_forecast",i) , forecast(get(paste0("mortality_model_fit",i)), 
                                                                 gc.order = get(paste0("gc.order_",i)),
                                                                 kt.method = "iarima",
                                                                 kt.order=unname(arimaorder(get(paste0("cv.arima.kt_",i)))),
                                                                 h = forecasting_horizon))
          
        }
        
        
        if (forecasting_horizon>1) {
          
          assign(paste0("mm_forecast_", i) , get(paste0("mortality_model_forecast",i))$rates[as.character(ages.fit), as.character(year.predict)])
          
        }else{
          
          assign(paste0("mm_forecast_", i) , get(paste0("mortality_model_forecast",i))$rates)
        }
        
        mse_2_ix <- mse_2_ix + (abs(get(paste0("muxt_actual_", i))-get(paste0("mm_forecast_", i))))/get(paste0("muxt_actual_", i))
        
      }
      
      mse_2 <- c(mse_2,sum(mse_2_ix))
      
    }
    
    
    
    out_model <- data.frame(
      forecasting_horizon = forecasting_horizon,
      model = c(model_option),
      mse_0 = c(sum(mse_0)/(length(ages.fit)*N_groups*length(prediction_horizon))),
      mse_1 = sum(mse_1)/(length(ages.fit)*N_groups*length(prediction_horizon)),
      mse_2 = sum(mse_2)/(length(ages.fit)*N_groups*length(prediction_horizon)))
    
    out <- rbind(out,
                 out_model)
    
  }}

fwrite(out,
       "C:\\Users\\gpitt\\Documents\\Postdoc\\Torino\\Mortality\\results\\nordic_validation.csv")
