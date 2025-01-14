# Fit mortality models ----

fit_mortality_models <- function(data_pp, years_fit, ages_fit) {
  l <- list()
  
  mortality_model_lc <- lc(link = "log")
  mortality_model_apc <- apc(link = "log")
  mortality_model_rh <- rh(link = "log")

  # Lee-Carter
  
  model_option <- "lc"
  
  assign("mortality_model", get(paste0("mortality_model_", model_option)))
  
  datahat <- data_pp$datahat
  list_of_extra_exposures <- data_pp$list_of_extra_exposures
  
  mortality_model_fit <- fit(
    mortality_model,
    data = datahat,
    years.fit = years_fit,
    ages.fit = ages_fit,
    verbose = FALSE,
    list_of_extra_exposures = list_of_extra_exposures
  )
  
  
  l[[model_option]] <- mortality_model_fit
  
  start.ax = mortality_model_fit$ax
  start.bx = mortality_model_fit$bx
  start.kt = mortality_model_fit$kt
  
  # Age-Period-Cohort
  
  model_option <- "apc"
  
  assign("mortality_model", get(paste0("mortality_model_", model_option)))
  
  mortality_model_fit <- fit(
    mortality_model,
    data = datahat,
    start.ax = start.ax,
    start.kt = start.kt,
    years.fit = years_fit,
    ages.fit = ages_fit,
    list_of_extra_exposures = list_of_extra_exposures,
    verbose = FALSE
  )
  
  
  l[[model_option]] <- mortality_model_fit
  
  # Renshaw- Habermann
  
  model_option <- "rh"
  
  assign("mortality_model", get(paste0("mortality_model_", model_option)))
  
  mortality_model_fit <- fit(
    mortality_model,
    data = datahat,
    start.ax = start.ax,
    start.kt = start.kt,
    years.fit = years_fit,
    ages.fit = ages_fit,
    list_of_extra_exposures = list_of_extra_exposures,
    verbose = FALSE
  )
  
  
  l[[model_option]] <- mortality_model_fit

  
  return(l)
  
  
}

# Fit credibility models ----

fit_and_predict_total_model <- function(data,
                                        data_pp,
                                        N_groups,
                                        mortality_models_fit,
                                        years_fit,
                                        ages_fit,
                                        scenario = 0,
                                        forecasting_horizon) {
  out <- list()
  observed_rates <- list()
  
  data_pp_2 <- data_preprocessing(
    data = data,
    N_groups = N_groups,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon)),
    scenario = scenario
  )
  
  for (model_option in names(mortality_models_fit)) {
    l <- list()
    
    mortality_model_fit <- mortality_models_fit[[model_option]]
    
    year.predict <- max(years_fit) + forecasting_horizon
    
    cv.arima.kt <- auto.arima(as.numeric(mortality_model_fit$kt), ic = "bic")
    
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
    
    muxt_hat <- mortality_model_forecast$fitted
    
    if (forecasting_horizon > 1) {
      muxt_hat_predicted <- mortality_model_forecast$rates[, as.character(last(years_fit) + forecasting_horizon)]
      
      
    } else{
      muxt_hat_predicted <- mortality_model_forecast$rates
      
      
    }
    
    tmp_occ <- data_pp$datahat$Dxt
    tmp_exp <- data_pp$datahat$Ext
    
    # C1 <- apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp * muxt_hat, 1, sum, na.rm =
    #                                                   T)
    # Fxt_1 <- tmp_occ / tmp_exp
    # tmpvar1 <- ((Fxt_1 - muxt_hat)^2) / (muxt_hat^2)
    # tmpvar1[tmpvar1 == 1] <- NA
    # varthetax_1 <- apply(tmpvar1, 1, mean, na.rm = T)
    # Z_1 <- 1 / (1 + varthetax_1 * apply(tmp_exp * muxt_hat, 1, sum, na.rm =
    #                                       T))
    #
    #
    subgroups_data <- data_pp$list_of_extra_exposures
    #
    #
    # l[['C1']] <- C1
    # l[['Fxt_1']] <- Fxt_1
    # l[['Z_1']] <- Z_1
    # l[['varthetax_1']] <- varthetax_1
    l[['muhat1']] <- muxt_hat_predicted
    observed_rates[['muxt_actual_1']] <- data_pp_2$datahat$Dxt[, (length(years_fit) + forecasting_horizon)] / data_pp_2$datahat$Ext[, (length(years_fit) + forecasting_horizon)]
    
    
    for (i in 2:N_groups) {
      tmp_occ <- subgroups_data[[i - 1]]$Dxt
      tmp_exp <- subgroups_data[[i - 1]]$Ext
      
      # assign(
      #   paste0("C", i),
      #   apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp *
      #                                               muxt_hat, 1, sum, na.rm = T)
      # )
      # assign(paste0("Fxt_", i), tmp_occ / tmp_exp)
      # assign(paste0("tmpvar", i), ((get(
      #   paste0("Fxt_", i)
      # ) - muxt_hat)^2) /
      #   (muxt_hat^2))
      #
      # tmp <- get(paste0("tmpvar", i))
      # tmp[tmp == 1] = NA
      # assign(paste0("tmpvar", i), tmp)
      # assign(paste0("varthetax_", i), apply(get(paste0("tmpvar", i)), 1, mean, na.rm =
      #                                         T))
      #
      # assign(paste0("Z_", i), 1 / (1 + get(paste0(
      #   "varthetax_", i
      # )) *
      #   apply(tmp_exp * muxt_hat, 1, sum, na.rm = T)))
      #
      # l[[paste0("C", i)]] <- get(paste0("C", i))
      # l[[paste0("Fxt_", i)]] <- get(paste0("Fxt_", i))
      # l[[paste0("Z_", i)]] <- get(paste0("Z_", i))
      # l[[paste0("varthetax_", i)]] <- get(paste0("varthetax_", i))
      
      l[[paste0("muhat", i)]] <-  muxt_hat_predicted
      
      subgroups_forecast <- data_pp_2$list_of_extra_exposures[[i - 1]]
      
      observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) + forecasting_horizon)]
      
      
    }
    
    observed_rates[['full_pp_data']] <- data_pp_2
    out[['actual_data']] <- observed_rates
    out[[model_option]] <- l
    
  }
  
  
  
  return(out)
}




fit_and_predict_credibility_models <- function(data,
                                               data_pp,
                                               N_groups,
                                               mortality_models_fit,
                                               years_fit,
                                               ages_fit,
                                               scenario = 0,
                                               forecasting_horizon) {
  out <- list()
  observed_rates <- list()
  
  data_pp_2 <- data_preprocessing(
    data = data,
    N_groups = N_groups,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon)),
    scenario = scenario
  )
  
  for (model_option in names(mortality_models_fit)) {
    l <- list()
    
    mortality_model_fit <- mortality_models_fit[[model_option]]
    
    year.predict <- max(years_fit) + forecasting_horizon
    
    cv.arima.kt <- auto.arima(as.numeric(mortality_model_fit$kt), ic = "bic")
    
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
    
    muxt_hat <- mortality_model_forecast$fitted
    
    if (forecasting_horizon > 1) {
      muxt_hat_predicted <- mortality_model_forecast$rates[, as.character(last(years_fit) + forecasting_horizon)]
      
      
    } else{
      muxt_hat_predicted <- mortality_model_forecast$rates
      
      
    }
    
    tmp_occ <- data_pp$datahat$Dxt
    tmp_exp <- data_pp$datahat$Ext
    
    C1 <- apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp * muxt_hat, 1, sum, na.rm =
                                                      T)
    Fxt_1 <- tmp_occ / tmp_exp
    tmpvar1 <- ((Fxt_1 - muxt_hat)^2) / (muxt_hat^2)
    tmpvar1[tmpvar1 == 1] <- NA
    varthetax_1 <- apply(tmpvar1, 1, mean, na.rm = T)
    Z_1 <- 1 / (1 + varthetax_1 * apply(tmp_exp * muxt_hat, 1, sum, na.rm =
                                          T))
    
    
    subgroups_data <- data_pp$list_of_extra_exposures
    
    
    l[['C1']] <- C1
    l[['Fxt_1']] <- Fxt_1
    l[['Z_1']] <- Z_1
    l[['varthetax_1']] <- varthetax_1
    l[['muhat1']] <- Z_1 * muxt_hat_predicted + (1 - Z_1) * C1 * muxt_hat_predicted
    
    # observed_rates[['muxt_actual_1']] <- data_pp_2$datahat$Dxt[, (length(years_fit) +
    #                                                                 1):(length(years_fit) + forecasting_horizon)] / data_pp_2$datahat$Ext[, (length(years_fit) +
    #                                                                                                                                            1):(length(years_fit) + forecasting_horizon)]
    observed_rates[['muxt_actual_1']] <- data_pp_2$datahat$Dxt[, (length(years_fit) + forecasting_horizon)] / data_pp_2$datahat$Ext[, (length(years_fit) + forecasting_horizon)]
    
    for (i in 2:N_groups) {
      tmp_occ <- subgroups_data[[i - 1]]$Dxt
      tmp_exp <- subgroups_data[[i - 1]]$Ext
      
      assign(
        paste0("C", i),
        apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp *
                                                    muxt_hat, 1, sum, na.rm = T)
      )
      assign(paste0("Fxt_", i), tmp_occ / tmp_exp)
      assign(paste0("tmpvar", i), ((get(
        paste0("Fxt_", i)
      ) - muxt_hat)^2) /
        (muxt_hat^2))
      
      tmp <- get(paste0("tmpvar", i))
      tmp[tmp == 1] = NA
      assign(paste0("tmpvar", i), tmp)
      assign(paste0("varthetax_", i), apply(get(paste0("tmpvar", i)), 1, mean, na.rm =
                                              T))
      
      assign(paste0("Z_", i), 1 / (1 + get(paste0(
        "varthetax_", i
      )) *
        apply(tmp_exp * muxt_hat, 1, sum, na.rm = T)))
      
      l[[paste0("C", i)]] <- get(paste0("C", i))
      l[[paste0("Fxt_", i)]] <- get(paste0("Fxt_", i))
      l[[paste0("Z_", i)]] <- get(paste0("Z_", i))
      l[[paste0("varthetax_", i)]] <- get(paste0("varthetax_", i))
      
      l[[paste0("muhat", i)]] <- get(paste0("Z_", i)) * muxt_hat_predicted + (1 - get(paste0("Z_", i))) * get(paste0("C", i)) * muxt_hat_predicted
      
      subgroups_forecast <- data_pp_2$list_of_extra_exposures[[i - 1]]
      
      # observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) +
      #                                                                            1):(length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) +
      #                                                                                                                                                        1):(length(years_fit) + forecasting_horizon)]
      #
      observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) + forecasting_horizon)]
      
    }
    
    observed_rates[['full_pp_data']] <- data_pp_2
    out[['actual_data']] <- observed_rates
    out[[model_option]] <- l
    
  }
  
  
  
  return(out)
}


fit_and_predict_separate_models <- function(data,
                                            data_pp,
                                            N_groups,
                                            mortality_models_fit,
                                            years_fit,
                                            ages_fit,
                                            scenario = 0,
                                            forecasting_horizon) {
  out <- list()
  observed_rates <- list()
  
  data_pp_2 <- data_preprocessing(
    data = data,
    N_groups = N_groups,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon)),
    scenario = scenario
  )
  
  for (model_option in names(mortality_models_fit)) {
    l <- list()
    i = 1
    popspecdata <- structure(
      list(
        Dxt = data_pp$datahat$Dxt[as.character(ages_fit), as.character(years_fit)],
        Ext = data_pp$datahat$Ext[as.character(ages_fit), as.character(years_fit)],
        ages = ages_fit,
        years = years_fit,
        type = 'central',
        series = "male",
        label =  paste0("superpop", i)
      ),
      class = "StMoMoData"
    )
    
    mortality_model_lc <- lc(link = "log")
    mortality_model_apc <- apc(link = "log")
    mortality_model_rh <- rh(link = "log")
    
    assign("mortality_model", get(paste0("mortality_model_", model_option)))
    
    assign(
      paste0("mortality_model_fit", i),
      fit(
        mortality_model,
        data = popspecdata,
        start.ax = mortality_models_fit[[model_option]]$ax,
        start.bx = mortality_models_fit[[model_option]]$bx,
        start.kt = mortality_models_fit[[model_option]]$kt,
        start.gc = mortality_models_fit[[model_option]]$gc,
        years.fit = years_fit,
        ages.fit = ages_fit,
        verbose = FALSE
      )
    )
    
    if(!get(paste0("mortality_model_fit", i))$conv){
      
      assign(
        paste0("mortality_model_fit", i),
        mortality_models_fit[[model_option]]
      )
      
      
    }
    
    assign(paste0("cv.arima.kt_", i),
           auto.arima(as.numeric(get(
             paste0("mortality_model_fit", i)
           )$kt), ic = "bic"))
    
    
    if (model_option != "lc") {
      assign(paste0("cv.arima.gc_", i),
             auto.arima(as.numeric(get(
               paste0("mortality_model_fit", i)
             )$gc), ic = "bic"))
      assign(paste0("gc.order_", i), unname(arimaorder(get(
        paste0("cv.arima.gc_", i)
      ))))
      
    } else{
      assign(paste0("gc.order_", i), c(1, 1, 0))
    }
    
    
    
    error_check <- tryCatch(
      assign(
        paste0("mortality_model_forecast", i) ,
        forecast(
          get(paste0("mortality_model_fit", i)),
          gc.order = get(paste0("gc.order_", i)),
          kt.method = "iarima",
          kt.order =
            unname(arimaorder(get(
              paste0("cv.arima.kt_", i)
            ))),
          h = forecasting_horizon
        )
      ),
      error = function(e) {
        return("mrw")
        
        
      }
    )
    
    
    if (is.character(error_check)) {
      assign(
        paste0("mortality_model_forecast", i) ,
        forecast(get(paste0(
          "mortality_model_fit", i
        )), h = forecasting_horizon)
      )
    } else{
      assign(
        paste0("mortality_model_forecast", i) ,
        forecast(
          get(paste0("mortality_model_fit", i)),
          gc.order = get(paste0("gc.order_", i)),
          kt.method = "iarima",
          kt.order =
            unname(arimaorder(get(
              paste0("cv.arima.kt_", i)
            ))),
          h = forecasting_horizon
        )
      )
      
    }
    
    
    if (forecasting_horizon > 1) {
      assign(paste0("muhat", i) , get(paste0("mortality_model_forecast", i))$rates[, as.character(last(years_fit) + forecasting_horizon)])
      
    } else{
      assign(paste0("muhat", i) , get(paste0("mortality_model_forecast", i))$rates)
    }
    
    l[['muhat1']] <- muhat1
    observed_rates[['muxt_actual_1']] <- data_pp_2$datahat$Dxt[, (length(years_fit) + forecasting_horizon)] / data_pp_2$datahat$Ext[, (length(years_fit) + forecasting_horizon)]
    
    
    for (i in 2:N_groups) {
      # browser()
      popspecdata <- structure(
        list(
          Dxt = data_pp$list_of_extra_exposures[[i - 1]]$Dxt[as.character(ages_fit), as.character(years_fit)],
          Ext = data_pp$list_of_extra_exposures[[i - 1]]$Ext[as.character(ages_fit), as.character(years_fit)],
          ages = ages_fit,
          years = years_fit,
          type = 'central',
          series = 'male',
          label =  paste0("superpop", i)
        ),
        class = "StMoMoData"
      )
      # browser()
      assign(
        paste0("mortality_model_fit", i),
        fit(
          mortality_model,
          data = popspecdata,
          start.ax = mortality_models_fit[[model_option]]$ax,
          start.bx = mortality_models_fit[[model_option]]$bx,
          start.kt = mortality_models_fit[[model_option]]$kt,
          start.gc = mortality_models_fit[[model_option]]$gc,
          years.fit = years_fit,
          ages.fit = ages_fit,
          verbose = FALSE
        )
      )
      
      if(!get(paste0("mortality_model_fit", i))$conv){
        
        assign(
          paste0("mortality_model_fit", i),
          fit(
            mortality_model,
            data = popspecdata,
            years.fit = years_fit,
            ages.fit = ages_fit,
            verbose = FALSE
          )
        )
        
        
      }
      
      if(!get(paste0("mortality_model_fit", i))$conv){
        
        assign(
          paste0("mortality_model_fit", i),
          mortality_models_fit[[model_option]]
        )
        
        
      }
      
      
      assign(paste0("cv.arima.kt_", i),
             auto.arima(as.numeric(get(
               paste0("mortality_model_fit", i)
             )$kt), ic = "bic"))
      
      
      if (model_option != "lc") {
        assign(paste0("cv.arima.gc_", i),
               auto.arima(as.numeric(get(
                 paste0("mortality_model_fit", i)
               )$gc), ic = "bic"))
        assign(paste0("gc.order_", i), unname(arimaorder(get(
          paste0("cv.arima.gc_", i)
        ))))
        
      } else{
        assign(paste0("gc.order_", i), c(1, 1, 0))
      }
      
      
      
      error_check <- tryCatch(
        assign(
          paste0("mortality_model_forecast", i) ,
          forecast(
            get(paste0("mortality_model_fit", i)),
            gc.order = get(paste0("gc.order_", i)),
            kt.method = "iarima",
            kt.order =
              unname(arimaorder(get(
                paste0("cv.arima.kt_", i)
              ))),
            h = forecasting_horizon
          )
        ),
        error = function(e) {
          return("mrw")
          
          
        }
      )
      
      
      if (is.character(error_check)) {
        assign(
          paste0("mortality_model_forecast", i) ,
          forecast(get(
            paste0("mortality_model_fit", i)
          ), h = forecasting_horizon)
        )
      } else{
        assign(
          paste0("mortality_model_forecast", i) ,
          forecast(
            get(paste0("mortality_model_fit", i)),
            gc.order = get(paste0("gc.order_", i)),
            kt.method = "iarima",
            kt.order =
              unname(arimaorder(get(
                paste0("cv.arima.kt_", i)
              ))),
            h = forecasting_horizon
          )
        )
        
      }
      
      
      if (forecasting_horizon > 1) {
        assign(paste0("muhat", i) , get(paste0("mortality_model_forecast", i))$rates[, as.character(last(years_fit) + forecasting_horizon)])
        
      } else{
        assign(paste0("muhat", i) , get(paste0("mortality_model_forecast", i))$rates)
      }
      
      
      l[[paste0("muhat", i)]] <- get(paste0("muhat", i))
      
      subgroups_forecast <- data_pp_2$list_of_extra_exposures[[i - 1]]
      
      observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) + forecasting_horizon)]
      
      
    }
    
    
    
    observed_rates[['full_pp_data']] <- data_pp_2
    out[['actual_data']] <- observed_rates
    out[[model_option]] <- l
    
  }
  
  
  
  return(out)
}



# Assess mortality model performance ----

poisson_nll <- function(occurrence, exposure, muxt_hat) {
  out <- -sum(occurrence * log(exposure * muxt_hat) - exposure * muxt_hat -
                lfactorial(occurrence))
  
  return(out)
  
}



model_assessment <- function(model_fit_and_prediction,
                             N_groups,
                             years_fit,
                             forecasting_horizon) {
  out <- list()
  
  relative_mse <- 0
  denominator <- 0
  poisson_nll_out <- 0
  tmp_actual_data <- model_fit_and_prediction[['actual_data']]
  
  
  for (model_ix in c("lc", "apc","rh")) {
    tmp_model <- model_fit_and_prediction[[model_ix]]
    
    l <- list()
    
    for (i in 1:N_groups) {
      tmp_group_string_actual <- paste0("muxt_actual_", i)
      tmp_group_string_fitted <- paste0("muhat", i)
      
      
      condition_selection <- !is.na(tmp_actual_data[[tmp_group_string_actual]]) &
        tmp_actual_data[[tmp_group_string_actual]] > 0
      
      denominator <- denominator + sum(condition_selection)
      
      # browser()
      relative_mse <- relative_mse + sum(((abs(
        tmp_model[[tmp_group_string_fitted]] - tmp_actual_data[[tmp_group_string_actual]]
      )) / tmp_actual_data[[tmp_group_string_actual]])[condition_selection])
      
      if (i < 2) {
        poisson_nll_out <- poisson_nll_out + poisson_nll(
          occurrence = tmp_actual_data$full_pp_data$datahat$Dxt[condition_selection, length(years_fit) +
                                                                  forecasting_horizon],
          exposure = tmp_actual_data$full_pp_data$datahat$Ext[condition_selection, length(years_fit) +
                                                                forecasting_horizon],
          muxt_hat = tmp_model[[tmp_group_string_fitted]][condition_selection]
        )
        
      } else{
        poisson_nll_out <- poisson_nll_out + poisson_nll(
          occurrence = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i - 1]][['Dxt']][condition_selection, length(years_fit) +
                                                                                                forecasting_horizon] ,
          exposure = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i -
                                                                             1]][['Ext']][condition_selection, length(years_fit) +
                                                                                            forecasting_horizon],
          muxt_hat = tmp_model[[tmp_group_string_fitted]][condition_selection]
        )
        
      }
      
      
      
    }
    
    l[['oos_deviance']] <- poisson_nll_out
    l[['error']] <- relative_mse / (denominator)
    
    out[[model_ix]] <- l
    
  }
  
  
  return(out)
  
  
  
}
