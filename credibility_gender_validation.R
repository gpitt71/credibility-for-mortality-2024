rm(list = ls())
library(StMoMo)
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(demography)

# Full model ----

G1.mx <- hmd.mx(country = "USA",
                "gabriele.pittarello@uniroma1.it",
                "pm869fqW8nxozS8!")

G2.mx <- hmd.mx(country = "USA",
                "gabriele.pittarello@uniroma1.it",
                "pm869fqW8nxozS8!")

ages.fit <- 50:89
years_fit_basic <- 1951:2003

series_code = 'male'
series_code_2 = 'female'

mortality_model_lc <- lc(link = "log")
mortality_model_apc <- apc(link = "log")
mortality_model_rh <- rh(link = "log")

N_groups = 2


out <- NULL 
for(forecasting_horizon in c(1,5,12)){
  
  prediction_horizon <- 0:(15 - forecasting_horizon)

  for (model_option in c("lc", "apc", "rh")) {
    

    mse_0 <- mse_1 <- mse_2 <- NULL

    assign("mortality_model", get(paste0("mortality_model_", model_option)))
    
    for (ix in prediction_horizon) {
      
      years.fit <- min(years_fit_basic):max(years_fit_basic + ix)
      
      E1 <- G1.mx$pop[[series_code]][as.character(ages.fit), as.character(years.fit)]
      E2 <- G2.mx$pop[[series_code_2]][as.character(ages.fit), as.character(years.fit)]
      
      D1 <- E1  * G1.mx$rate[[series_code]][as.character(ages.fit), as.character(years.fit)]
      D2 <- E2  * G2.mx$rate[[series_code_2]][as.character(ages.fit), as.character(years.fit)]
      
      Ext <- E1 #+ E2
      Dxt <-  D1 #+ D2
      
      
      datahat <- structure(
        list(
          Dxt = Dxt,
          Ext = Ext,
          ages = ages.fit,
          years = years.fit,
          type = 'central',
          series = series_code,
          label = 'total'
        ),
        class = "StMoMoData"
      )
      
      
      list_of_extra_exposures <- list()
      
      list_of_extra_exposures[[1]] <- list(Dxt=D2,
                                           Ext=E2)
      
      mortality_model_fit <- fit(mortality_model,
                        data = datahat,
                        years.fit = years.fit,
                        ages.fit = ages.fit,
                        list_of_extra_exposures=list_of_extra_exposures)
      
      year.predict <- max(years.fit) + forecasting_horizon
      

      cv.arima.kt <- auto.arima(as.numeric(mortality_model_fit$kt),
                                ic="bic")
      
      if(model_option != "lc"){
        cv.arima.gc <- auto.arima(as.numeric(mortality_model_fit$gc),
                                  ic="bic")
        gc.order <- unname(arimaorder(cv.arima.gc))
      }else{
        
        gc.order <- c(1,1,0)
      }
      

      mortality_model_forecast <- forecast(mortality_model_fit, 
                                           kt.method = "iarima",
                                           gc.order=gc.order,
                                           kt.order=unname(arimaorder(cv.arima.kt)),
                                           h = forecasting_horizon)
      
      
      

      muxt_actual_1 <- G1.mx$rate[[series_code]][as.character(ages.fit), as.character(year.predict)]
      muxt_actual_2 <- G2.mx$rate[[series_code_2]][as.character(ages.fit), as.character(year.predict)]
      
      

      
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
      
      C1 <- apply(D1, 1, sum) / apply(E1 * muxt_hat, 1, sum)
      C2 <- apply(D2, 1, sum) / apply(E2 * muxt_hat, 1, sum)

      
      Fxt_1 <- D1 / E1
      Fxt_2 <- D2 / E2
      
      varthetax_1 <- apply(((Fxt_1 - muxt_hat) ^ 2) / (muxt_hat ^ 2), 1, mean)
      varthetax_2 <- apply(((Fxt_2 - muxt_hat) ^ 2) / (muxt_hat ^ 2), 1, mean)
      
      Z_1 <- 1 / (1 + varthetax_1 * apply(E1 * muxt_hat, 1, sum))

      Z_2 <- 1 / (1 + varthetax_2 * apply(E2 * muxt_hat, 1, sum))

      muhat1 <- Z_1 * muxt_hat_predicted + (1 - Z_1) * C1 * muxt_hat_predicted
      muhat2 <- Z_2 * muxt_hat_predicted + (1 - Z_2) * C2 * muxt_hat_predicted
      

      mse_1_ix <- 0

      for (i in 1:N_groups) {
        
        mse_1_ix <- mse_1_ix + (abs(get(paste0("muxt_actual_", i))-get(paste0("muhat", i))))#/get(paste0("muxt_actual_", i))

      }
      
      mse_1 <- c(mse_1,sum(mse_1_ix))

      # Separated model ----
      
      G1data <- StMoMoData(G1.mx, series = series_code, type = "central")
      mortality_model_fit1 <- fit(mortality_model,
                         data = G1data,
                         years.fit = years.fit,
                         ages.fit = ages.fit)
      
      cv.arima.kt_1 <- auto.arima(as.numeric(mortality_model_fit1$kt),
                                  ic="bic")
      
      G2data <- StMoMoData(G2.mx, series = series_code_2, type = "central")
      mortality_model_fit2 <- fit(mortality_model,
                         data = G2data,
                         years.fit = years.fit,
                         ages.fit = ages.fit)
      
      cv.arima.kt_2 <- auto.arima(as.numeric(mortality_model_fit2$kt),
                                  ic="bic")
      
      if(model_option != "lc"){
        cv.arima.gc_1 <- auto.arima(as.numeric(mortality_model_fit1$gc),
                                    ic="bic")
        gc.order_1 <- unname(arimaorder(cv.arima.gc_1))
        cv.arima.gc_2 <- auto.arima(as.numeric(mortality_model_fit2$gc),
                                    ic="bic")
        gc.order_2 <- unname(arimaorder(cv.arima.gc_2))

      }else{
        
        gc.order_1 <- gc.order_2 <- c(1,1,0)
      }
      
      
      
      mortality_model_forecast1 <- forecast(mortality_model_fit1, 
                                            gc.order = gc.order_1,
                                            kt.method = "iarima",
                                            kt.order=unname(arimaorder(cv.arima.kt_1)),
                                            h = forecasting_horizon)
      
      mortality_model_forecast2 <- forecast(mortality_model_fit2, 
                                            gc.order = gc.order_2,
                                            kt.method = "iarima",
                                            kt.order=unname(arimaorder(cv.arima.kt_2)),
                                            h = forecasting_horizon)
      
      
      if (forecasting_horizon>1) {
        mm_forecast_1 <- mortality_model_forecast1$rates[as.character(ages.fit), as.character(year.predict)]

        mm_forecast_2 <- mortality_model_forecast2$rates[as.character(ages.fit), as.character(year.predict)]

        }else{
        mm_forecast_1 <- mortality_model_forecast1$rates
        mm_forecast_2 <- mortality_model_forecast2$rates
      }
      
      mse_2_ix <- 0

      for(i in 1:N_groups){
        mse_2_ix <- mse_2_ix + (abs(get(paste0("muxt_actual_", i))-get(paste0("mm_forecast_", i))))#/get(paste0("muxt_actual_", i))     

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
    
  }
}


fwrite(out,
       "C:\\Users\\gpitt\\Documents\\Postdoc\\Torino\\Mortality\\results\\gender_validation_l2.csv")
