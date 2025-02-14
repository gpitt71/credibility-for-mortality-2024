rm(list = ls())
library(StMoMo)
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(demography)

# Full model ----

DK.mx <- hmd.mx(country = "DNK",
                "gabriele.pittarello@uniroma1.it",
                "pm869fqW8nxozS8!")

SWE.mx <- hmd.mx(country = "SWE",
                 "gabriele.pittarello@uniroma1.it",
                 "pm869fqW8nxozS8!")

NOR.mx <- hmd.mx(country = "NOR",
                 "gabriele.pittarello@uniroma1.it",
                 "pm869fqW8nxozS8!")


ages.fit <- 50:89
years_fit_basic <- 1951:2003

series_code = 'male'

mortality_model_lc <- lc(link = "log")
mortality_model_apc <- apc(link = "log")
mortality_model_rh <- rh(link = "log")

out <- NULL

N_groups = 3

for (forecasting_horizon in c(1, 5, 12)) {
  prediction_horizon <- 0:(15 - forecasting_horizon)
  
  for (model_option in c("lc", "apc", "rh")) {
    mse_0 <- mse_1 <- mse_2 <- NULL
    
    assign("mortality_model", get(paste0("mortality_model_", model_option)))
    
    for (ix in prediction_horizon) {
      years.fit <- min(years_fit_basic):max(years_fit_basic + ix)
      
      E1 <- DK.mx$pop[[series_code]][as.character(ages.fit), as.character(years.fit)]
      E2 <- SWE.mx$pop[[series_code]][as.character(ages.fit), as.character(years.fit)]
      E3 <- NOR.mx$pop[[series_code]][as.character(ages.fit), as.character(years.fit)]
      
      D1 <- E1  * DK.mx$rate[[series_code]][as.character(ages.fit), as.character(years.fit)]
      D2 <- E2  * SWE.mx$rate[[series_code]][as.character(ages.fit), as.character(years.fit)]
      D3 <- E3  * NOR.mx$rate[[series_code]][as.character(ages.fit), as.character(years.fit)]
      
      Ext <- E1 #+ E2 + E3
      Dxt <-  D1 #+ D2 + D3
      
      list_of_extra_exposures <- list()
      
      for (i in 2:N_groups) {
        tmp_list <- list()
        
        tmp_list[['Dxt']] <- get(paste0("D", i))
        tmp_list[['Ext']] <- get(paste0("E", i))
        list_of_extra_exposures[[i - 1]] <- tmp_list
        
      }
      
      
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
      
      
      
      mortality_model_fit <- fit(
        mortality_model,
        data = datahat,
        years.fit = years.fit,
        ages.fit = ages.fit,
        list_of_extra_exposures = list_of_extra_exposures
      )
      
      year.predict <- max(years.fit) + forecasting_horizon
      
      cv.arima.kt <- auto.arima(as.numeric(mortality_model_fit$kt), ic =
                                  "bic")
      
      if (model_option != "lc") {
        cv.arima.gc <- auto.arima(as.numeric(mortality_model_fit$gc), ic = "bic")
        gc.order <- unname(arimaorder(cv.arima.gc))
      } else{
        gc.order <- c(1, 1, 0)
      }
      
      mortality_model_forecast <- forecast(
        mortality_model_fit,
        kt.method = "iarima",
        gc.order = gc.order,
        kt.order = unname(arimaorder(cv.arima.kt)),
        h = forecasting_horizon
      )
      
      
      muxt_actual_1 <- DK.mx$rate[[series_code]][as.character(ages.fit), as.character(year.predict)]
      muxt_actual_2 <- SWE.mx$rate[[series_code]][as.character(ages.fit), as.character(year.predict)]
      muxt_actual_3 <- NOR.mx$rate[[series_code]][as.character(ages.fit), as.character(year.predict)]
      
      if (forecasting_horizon > 1) {
        muxt_hat_predicted <- mortality_model_forecast$rates[as.character(ages.fit), as.character(year.predict)]
      } else{
        muxt_hat_predicted <- mortality_model_forecast$rates
      }
      
      
      mse_0_ix <- 0
      
      for (i in 1:N_groups) {
        mse_0_ix <- mse_0_ix + (abs((
          get(paste0("muxt_actual_", i)) - muxt_hat_predicted
        )))/get(paste0("muxt_actual_", i))
        
      }
      
      mse_0 <- c(mse_0, sum(mse_0_ix))
      
      # Credibility model -----
      
      muxt_hat <- mortality_model_forecast$fitted
      
      C1 <- apply(D1, 1, sum) / apply(E1 * muxt_hat, 1, sum)
      
      C2 <- apply(D2, 1, sum) / apply(E2 * muxt_hat, 1, sum)
      
      C3 <- apply(D3, 1, sum) / apply(E3 * muxt_hat, 1, sum)
      
      Fxt_1 <- D1 / E1
      Fxt_2 <- D2 / E2
      Fxt_3 <- D3 / E3
      
      varthetax_1 <- apply(((Fxt_1 - muxt_hat)^2) / (muxt_hat^2), 1, mean)
      varthetax_2 <- apply(((Fxt_2 - muxt_hat)^2) / (muxt_hat^2), 1, mean)
      varthetax_3 <- apply(((Fxt_3 - muxt_hat)^2) / (muxt_hat^2), 1, mean)
      
      Z_1 <- 1 / (1 + varthetax_1 * apply(E1 * muxt_hat, 1, sum))
      
      Z_2 <- 1 / (1 + varthetax_2 * apply(E2 * muxt_hat, 1, sum))
      
      Z_3 <- 1 / (1 + varthetax_3 * apply(E3 * muxt_hat, 1, sum))
      
      muhat1 <- Z_1 * muxt_hat_predicted + (1 - Z_1) * C1 * muxt_hat_predicted
      muhat2 <- Z_2 * muxt_hat_predicted + (1 - Z_2) * C2 * muxt_hat_predicted
      muhat3 <- Z_3 * muxt_hat_predicted + (1 - Z_3) * C3 * muxt_hat_predicted
      
      mse_1_ix <- 0
      
      for (i in 1:N_groups) {
        mse_1_ix <- mse_1_ix + (abs(get(paste0(
          "muxt_actual_", i
        )) - get(paste0("muhat", i))))/get(paste0("muxt_actual_", i))
        
      }
      
      mse_1 <- c(mse_1, sum(mse_1_ix))
      
      # Separated model ----
      
      DKdata <- StMoMoData(DK.mx, series = series_code, type = "central")
      mortality_model_fit1 <- fit(
        mortality_model,
        data = DKdata,
        years.fit = years.fit,
        ages.fit = ages.fit
      )
      
      cv.arima.kt_1 <- auto.arima(as.numeric(mortality_model_fit1$kt), ic =
                                    "bic")
      
      SWEdata <- StMoMoData(SWE.mx, series = series_code, type = "central")
      mortality_model_fit2 <- fit(
        mortality_model,
        data = SWEdata,
        years.fit = years.fit,
        ages.fit = ages.fit
      )
      cv.arima.kt_2 <- auto.arima(as.numeric(mortality_model_fit2$kt), ic =
                                    "bic")
      
      
      NORdata <- StMoMoData(NOR.mx, series = series_code, type = "central")
      
      mortality_model_fit3 <- fit(
        mortality_model,
        data = NORdata,
        years.fit = years.fit,
        ages.fit = ages.fit
      )
      cv.arima.kt_3 <- auto.arima(as.numeric(mortality_model_fit3$kt), ic =
                                    "bic")
      
      if (model_option != "lc") {
        cv.arima.gc_1 <- auto.arima(as.numeric(mortality_model_fit1$gc), ic = "bic")
        gc.order_1 <- unname(arimaorder(cv.arima.gc_1))
        cv.arima.gc_2 <- auto.arima(as.numeric(mortality_model_fit2$gc), ic =
                                      "bic")
        gc.order_2 <- unname(arimaorder(cv.arima.gc_2))
        cv.arima.gc_3 <- auto.arima(as.numeric(mortality_model_fit3$gc), ic =
                                      "bic")
        gc.order_3 <- unname(arimaorder(cv.arima.gc_3))
      } else{
        gc.order_1 <- gc.order_2 <- gc.order_3 <- c(1, 1, 0)
      }
      
      mortality_model_forecast1 <- forecast(
        mortality_model_fit1,
        gc.order = gc.order_1,
        kt.method = "iarima",
        kt.order = unname(arimaorder(cv.arima.kt_1)),
        h = forecasting_horizon
      )
      mortality_model_forecast2 <- forecast(
        mortality_model_fit2,
        gc.order = gc.order_2,
        kt.method = "iarima",
        kt.order = unname(arimaorder(cv.arima.kt_2)),
        h = forecasting_horizon
      )
      mortality_model_forecast3 <- forecast(
        mortality_model_fit3,
        gc.order = gc.order_3,
        kt.method = "iarima",
        kt.order = unname(arimaorder(cv.arima.kt_3)),
        ,
        h = forecasting_horizon
      )
      
      if (forecasting_horizon > 1) {
        mm_forecast_1 <- mortality_model_forecast1$rates[as.character(ages.fit), as.character(year.predict)]
        mm_forecast_2 <- mortality_model_forecast2$rates[as.character(ages.fit), as.character(year.predict)]
        mm_forecast_3 <- mortality_model_forecast3$rates[as.character(ages.fit), as.character(year.predict)]
      } else{
        mm_forecast_1 <- mortality_model_forecast1$rates
        mm_forecast_2 <- mortality_model_forecast2$rates
        mm_forecast_3 <- mortality_model_forecast3$rates
      }
      
      mse_2_ix <- 0
      
      for (i in 1:N_groups) {
        mse_2_ix <- mse_2_ix + (abs(get(paste0(
          "muxt_actual_", i
        )) - get(paste0(
          "mm_forecast_", i
        ))))/get(paste0("muxt_actual_", i))
        
      }
      
      mse_2 <- c(mse_2, sum(mse_2_ix))
      
    }
    
    
    
    out_model <- data.frame(
      forecasting_horizon = forecasting_horizon,
      model = c(model_option),
      mse_0 = c(sum(mse_0) / (
        length(ages.fit) * N_groups * length(prediction_horizon)
      )),
      mse_1 = sum(mse_1) / (length(ages.fit) * N_groups * length(prediction_horizon)),
      mse_2 = sum(mse_2) / (length(ages.fit) * N_groups * length(prediction_horizon))
    )
    
    out <- rbind(out, out_model)
    
  }
}

data.table::fwrite(
  out,
  paste0("~/GitHub/credibility-for-mortality-2024/output/scandinavian_rolling_window_",series_code,".csv")
)
