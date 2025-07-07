# Fit mortality models ----

fit_mortality_models <- function(data_pp,
                                 years_fit,
                                 ages_fit,
                                 separate_exposures=FALSE) {
  l <- list()

  mortality_model_lc <- lc(link = "log")
  mortality_model_apc <- apc(link = "log")
  mortality_model_rh <- rh(link = "log")

  # Lee-Carter

  model_option <- "lc"

  assign("mortality_model", get(paste0("mortality_model_", model_option)))

  datahat <- data_pp$datahat

  list_of_extra_exposures <- data_pp$list_of_extra_exposures

  if(separate_exposures){

  mortality_model_fit <- fit(
    mortality_model,
    data = datahat,
    years.fit = years_fit,
    ages.fit = ages_fit,
    verbose = FALSE
    ,
    list_of_extra_exposures = (list_of_extra_exposures)
  )
  }else{

    datahat <- data_pp$datahat_tot

    mortality_model_fit <- fit(
      mortality_model,
      data = datahat,
      years.fit = years_fit,
      ages.fit = ages_fit,
      verbose = FALSE
    )
  }


  l[[model_option]] <- mortality_model_fit

  start.ax = mortality_model_fit$ax
  start.bx = mortality_model_fit$bx
  start.kt = mortality_model_fit$kt

  # Age-Period-Cohort

  model_option <- "apc"

  assign("mortality_model", get(paste0("mortality_model_", model_option)))

  if(separate_exposures){

    mortality_model_fit <- fit(
      mortality_model,
      data = datahat,
      years.fit = years_fit,
      ages.fit = ages_fit,
      verbose = FALSE
      ,
      list_of_extra_exposures = (list_of_extra_exposures)
    )
  }else{

    datahat <- data_pp$datahat_tot

    mortality_model_fit <- fit(
      mortality_model,
      data = datahat,
      years.fit = years_fit,
      ages.fit = ages_fit,
      verbose = FALSE
    )
  }




  l[[model_option]] <- mortality_model_fit

  # Renshaw- Habermann

  model_option <- "rh"

  assign("mortality_model", get(paste0("mortality_model_", model_option)))

  if(separate_exposures){

    mortality_model_fit <- fit(
      mortality_model,
      data = datahat,
      years.fit = years_fit,
      ages.fit = ages_fit,
      verbose = FALSE
      ,
      list_of_extra_exposures = (list_of_extra_exposures)
    )
  }else{

    datahat <- data_pp$datahat_tot

    mortality_model_fit <- fit(
      mortality_model,
      data = datahat,
      years.fit = years_fit,
      ages.fit = ages_fit,
      verbose = FALSE
    )
  }

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
                                        bias=70,
                                        scenario = 2,
                                        forecasting_horizon) {
  out <- list()
  observed_rates <- list()

  data_pp_2 <- data_preprocessing_scenario_2(
    data = data,
    N_groups = N_groups,
    bias = bias,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon))
    # scenario = scenario
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

    l[['muhatsuperpop']] <- muxt_hat_predicted
    observed_rates[['muxt_actual_superpop']] <- data_pp_2$datahat_tot$Dxt[, (length(years_fit) + forecasting_horizon)] / data_pp_2$datahat_tot$Ext[, (length(years_fit) + forecasting_horizon)]


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
                                               bias=70,
                                               scenario = 2,
                                               smoothing_years=5,
                                               dependent_variance=TRUE,
                                               forecasting_horizon) {

  out <- list()
  observed_rates <- list()

  data_pp_2 <- data_preprocessing_scenario_2(
    data = data,
    N_groups = N_groups,
    bias = bias,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon))
    # scenario = scenario
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

    # browser()

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

    tmp_occ <- data_pp$datahat$Dxt[as.character(ages_fit), as.character(years_fit)]
    tmp_exp <- data_pp$datahat$Ext[as.character(ages_fit), as.character(years_fit)]

    C1 <- apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp * muxt_hat, 1, sum, na.rm =
                                                      T)
    # browser()

    fit1 <- rpart(C1~ages_fit)
    C1 <- predict(fit1, data.frame(ages_fit))


    # C1 <- k_years_smoothing(C1, k= smoothing_years)

    Fxt_1 <- tmp_occ / tmp_exp


    # tmpvar1[tmpvar1 == 1] <- NA
    # browser()

    if(dependent_variance){

      tmpvar1 <- apply(Fxt_1,1,sum,na.rm=T)
      tmpvar2 <- apply(muxt_hat,1,sum,na.rm=T)

      numerator <- (tmpvar1 - tmpvar2)^2 - apply(muxt_hat / tmp_exp, 1, sum, na.rm = T)

      denominator <- apply(muxt_hat,1,sum,na.rm=TRUE)^2

      varthetax_1 <- pmax(0,
                          numerator / denominator,
                          na.rm = T)

    }else{

      tmpvar1 <- ((Fxt_1 - muxt_hat)^2 - (muxt_hat / tmp_exp))

      varthetax_1 <- pmax(0,
                          apply(tmpvar1, 1, sum, na.rm = T) / apply(muxt_hat^2, 1, sum, na.rm = T),
                          na.rm = T)

}
    # varthetax_1 <- k_years_smoothing(varthetax_1, k= smoothing_years)

    fit1 <- rpart(varthetax_1~ages_fit)
    varthetax_1 <- predict(fit1, data.frame(ages_fit))

    Z_1 <- 1 / (1 + varthetax_1 * apply(tmp_exp * muxt_hat, 1, sum, na.rm =
                                          T))

    subgroups_data <- data_pp$list_of_extra_exposures


    l[['C1']] <- C1
    l[['Fxt_1']] <- Fxt_1
    l[['Z_1']] <- Z_1
    l[['varthetax_1']] <- varthetax_1
    l[['muhat1']] <- Z_1 * muxt_hat_predicted + (1 - Z_1) * C1 * muxt_hat_predicted
    l[['muhat1_MLE']] <- C1 * muxt_hat_predicted


    # term1 <- 2*(muxt_hat_predicted^2)*(varthetax_1)+muxt_hat_predicted^2
    #
    # term2 <- varthetax_1*(apply((tmp_exp^2) *( muxt_hat)^2, 1, sum, na.rm =
    #                   T)/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
    #                     T)^2)+(1/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
    #                                     T))
    #
    # l[['msep_1']] <- term1+term2


    # observed_rates[['muxt_actual_1']] <- data_pp_2$datahat$Dxt[, (length(years_fit) +
    #                                                                 1):(length(years_fit) + forecasting_horizon)] / data_pp_2$datahat$Ext[, (length(years_fit) +
    #                                                                                                                                            1):(length(years_fit) + forecasting_horizon)]
    observed_rates[['muxt_actual_1']] <- data_pp_2$datahat$Dxt[, (length(years_fit) + forecasting_horizon)] / data_pp_2$datahat$Ext[, (length(years_fit) + forecasting_horizon)]

    for (i in 2:N_groups) {
      tmp_occ <- subgroups_data[[i - 1]]$Dxt[as.character(ages_fit), as.character(years_fit)]
      tmp_exp <- subgroups_data[[i - 1]]$Ext[as.character(ages_fit), as.character(years_fit)]

      assign(
        paste0("C", i),
        apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp *
                                                    muxt_hat, 1, sum, na.rm = T)
      )

      fit1 <- rpart(get(paste0("C", i))~ages_fit)
      assign(paste0("C", i), predict(fit1, data.frame(ages_fit)))

      assign(paste0("Fxt_", i), tmp_occ / tmp_exp)


      if(dependent_variance){
        # browser()
        # tmpvar1 <- apply(Fxt_1,1,sum,na.rm=T)
        # tmpvar2 <- apply(muxt_hat,1,sum,na.rm=T)

        assign(paste0("tmpvar", i), #((

          (apply(get(
          paste0("Fxt_", i)
        ), 1, sum, na.rm = T) - apply(muxt_hat,1,sum,na.rm=T))^2 - apply(muxt_hat / tmp_exp, 1, sum, na.rm = T)
        )
        #- muxt_hat)^2 - (muxt_hat / tmp_exp)))

        # numerator <- (tmpvar1 - tmpvar2)^2 - apply(muxt_hat / tmp_exp, 1, sum, na.rm = T)

        denominator <- apply(muxt_hat,1,sum,na.rm=TRUE)^2

        assign(paste0("varthetax_", i), pmax(
          0,
          get(paste0("tmpvar", i))/ denominator,
          na.rm = T
        ))

      }else{


        assign(paste0("tmpvar", i), ((get(
          paste0("Fxt_", i)
        ) - muxt_hat)^2 - (muxt_hat / tmp_exp)))

        tmp <- get(paste0("tmpvar", i))
        # tmp[tmp == 1] = NA


      assign(paste0("tmpvar", i), tmp)
      assign(paste0("varthetax_", i), pmax(
        0,
        apply(get(paste0("tmpvar", i)), 1, mean, na.rm =
                T) / apply(muxt_hat^2, 1, mean, na.rm =
                             T),
        na.rm = T
      ))

      }

      fit1 <- rpart(get(paste0("varthetax_", i))~ages_fit)
      assign(paste0("varthetax_", i), predict(fit1, data.frame(ages_fit)))

      assign(paste0("Z_", i), 1 / (1 + get(paste0(
        "varthetax_", i
      )) *
        apply(tmp_exp * muxt_hat, 1, sum, na.rm = T)))




      l[[paste0("C", i)]] <- get(paste0("C", i))
      l[[paste0("Fxt_", i)]] <- get(paste0("Fxt_", i))
      l[[paste0("Z_", i)]] <- get(paste0("Z_", i))
      l[[paste0("varthetax_", i)]] <- get(paste0("varthetax_", i))

      l[[paste0("muhat", i)]] <- get(paste0("Z_", i)) * muxt_hat_predicted + (1 - get(paste0("Z_", i))) * get(paste0("C", i)) * muxt_hat_predicted
      l[[paste0("muhat", i,"_MLE")]] <- get(paste0("C", i)) * muxt_hat_predicted


      # term1 <- 2*(muxt_hat_predicted^2)*(get(paste0("varthetax_", i)))+muxt_hat_predicted^2
      #
      # term2 <- get(paste0("varthetax_", i))*(apply((tmp_exp^2) *( muxt_hat)^2, 1, sum, na.rm =
      #                               T)/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
      #                                           T)^2)+(1/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
      #                                                             T))
      #
      # l[[paste0('msep_',i)]] <- term1+term2

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
                                            bias=70,
                                            scenario = 2,
                                            forecasting_horizon) {

  out <- list()
  observed_rates <- list()
  failed_fits <- data.frame(group=9999,
                            failed_fit=99999,
                            model="bla")



  data_pp_2 <- data_preprocessing_scenario_2(
    data = data,
    N_groups = N_groups,
    bias = bias,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon))
    # scenario = scenario
  )

  for (model_option in names(mortality_models_fit)) {
    tmp_failed=FALSE
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

    error_check <- tryCatch(
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
      ),
      error = function(e) {
        return("failed")


      }
    )

    if (is.character(error_check)) {
      assign(paste0("mortality_model_fit", i),
             mortality_models_fit[[model_option]])

      tmp_failed=TRUE



    } else{
      assign(paste0("mortality_model_fit", i), error_check)


      if (!get(paste0("mortality_model_fit", i))$conv) {
        assign(paste0("mortality_model_fit", i),
               mortality_models_fit[[model_option]])

        tmp_failed=TRUE



      }



    }

    if(tmp_failed){failed_fits = rbind(failed_fits, data.frame(group=i,failed_fit=1,model=model_option))}

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


    l[[paste0("model", i)]] <- get(paste0("mortality_model_fit", i))

    for (i in 2:N_groups) {

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






      error_check <- tryCatch(
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
        ),
        error = function(e) {
          return("failed")


        }
      )

      tmp_failed=FALSE

      if (is.character(error_check)) {

        # browser()

        assign(paste0("mortality_model_fit", i),
               mortality_models_fit[[model_option]])

        tmp_failed=TRUE



      } else{

        # browser()

        assign(paste0("mortality_model_fit", i), error_check)


        if (!get(paste0("mortality_model_fit", i))$conv) {
          assign(paste0("mortality_model_fit", i),
                 mortality_models_fit[[model_option]])

          tmp_failed=TRUE


        }

      }





      if (!get(paste0("mortality_model_fit", i))$conv) {
        assign(paste0("mortality_model_fit", i),
               mortality_models_fit[[model_option]])


        tmp_failed=TRUE
      }

      if(tmp_failed){failed_fits = rbind(failed_fits, data.frame(group=i,failed_fit=1,model=model_option))}


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

      l[[paste0("model", i)]] <- get(paste0("mortality_model_fit", i))

      subgroups_forecast <- data_pp_2$list_of_extra_exposures[[i - 1]]

      observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) + forecasting_horizon)]


    }



    observed_rates[['full_pp_data']] <- data_pp_2
    out[['actual_data']] <- observed_rates
    out[[model_option]] <- l

  }


  if(nrow(failed_fits)>1){

  output_folder <- "~/credibility-for-mortality-2024/output/"


  timestamp <- format(Sys.time(), "%H_%M_%S_%d")
  file_name <- paste0("failed_fit_", timestamp, ".txt")


  file_path <- file.path(output_folder, file_name)


  fwrite(failed_fits[-1,], file = file_path, sep = "\t")
  }



  return(out)
}



# Assess mortality model performance ----

poisson_nll <- function(occurrence, exposure, muxt_hat) {
  out <- -sum(occurrence * log(exposure * muxt_hat) - exposure * muxt_hat -
                lfactorial(occurrence))

  return(out)

}


# Model assessment ----

model_assessment <- function(model_fit_and_prediction,
                             N_groups,
                             years_fit,
                             forecasting_horizon,
                             bias=15,
                             subset_rates=NULL) {
  out <- list()

  relative_mse <- 0
  denominator <- 0
  poisson_dev_out <- 0
  tmp_actual_data <- model_fit_and_prediction[['actual_data']]


  for (model_ix in c("lc", "apc", "rh")) {
    tmp_model <- model_fit_and_prediction[[model_ix]]

    l <- list()


    for (i in 1:N_groups) {


      tmp_group_string_actual <- paste0("muxt_actual_", i)
      tmp_group_string_fitted <- paste0("muhat", i)

      if(!is.null(subset_rates)){

        subsetting_condition <- rep(FALSE, length(tmp_actual_data[[tmp_group_string_actual]]))

        subsetting_condition[subset_rates-bias+1] <- TRUE

      }else{

        subsetting_condition <- rep(TRUE, length(tmp_actual_data[[tmp_group_string_actual]]))

      }



      # browser()
      condition_selection <- !is.na(tmp_actual_data[[tmp_group_string_actual]]) &
        tmp_actual_data[[tmp_group_string_actual]] > 0 &subsetting_condition

      denominator <- denominator + sum(condition_selection)


      relative_mse <- relative_mse + sum(((abs(
        tmp_model[[tmp_group_string_fitted]] - tmp_actual_data[[tmp_group_string_actual]]
      )) / tmp_actual_data[[tmp_group_string_actual]])[condition_selection])

      relative_mse_g <- (sum(((abs(
        tmp_model[[tmp_group_string_fitted]] - tmp_actual_data[[tmp_group_string_actual]]
      )) / tmp_actual_data[[tmp_group_string_actual]])[condition_selection]))/sum(condition_selection)

      if (i < 2) {

        models_ll <- -poisson_nll(
          occurrence = tmp_actual_data$full_pp_data$datahat$Dxt[condition_selection, length(years_fit) +
                                                                  forecasting_horizon],
          exposure = tmp_actual_data$full_pp_data$datahat$Ext[condition_selection, length(years_fit) +
                                                                forecasting_horizon],
          muxt_hat = tmp_model[[tmp_group_string_fitted]][condition_selection]
        )

        saturated_model_ll <- -poisson_nll(
          occurrence = tmp_actual_data$full_pp_data$datahat$Dxt[condition_selection, length(years_fit) +
                                                                  forecasting_horizon],
          exposure = tmp_actual_data$full_pp_data$datahat$Ext[condition_selection, length(years_fit) +
                                                                forecasting_horizon],
          muxt_hat = tmp_actual_data$full_pp_data$datahat$Dxt[condition_selection, length(years_fit) +
                                                                forecasting_horizon]/tmp_actual_data$full_pp_data$datahat$Ext[condition_selection, length(years_fit) +
                                                                                                                                forecasting_horizon]
        )

        poisson_dev_out <- poisson_dev_out  -2*(models_ll-saturated_model_ll)

      } else{


        models_ll <- -poisson_nll(
          occurrence = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i - 1]][['Dxt']][condition_selection, length(years_fit) +
                                                                                                forecasting_horizon] ,
          exposure = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i -
                                                                             1]][['Ext']][condition_selection, length(years_fit) +
                                                                                            forecasting_horizon],
          muxt_hat = tmp_model[[tmp_group_string_fitted]][condition_selection]
        )

        saturated_model_ll <- -poisson_nll(
          occurrence = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i - 1]][['Dxt']][condition_selection, length(years_fit) +
                                                                                                forecasting_horizon] ,
          exposure = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i -
                                                                             1]][['Ext']][condition_selection, length(years_fit) +
                                                                                            forecasting_horizon],
          muxt_hat = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i - 1]][['Dxt']][condition_selection, length(years_fit) +
                                                                                              forecasting_horizon]/tmp_actual_data$full_pp_data$list_of_extra_exposures[[i -
                                                                                                                                                                           1]][['Ext']][condition_selection, length(years_fit) +
                                                                                                                                                                                          forecasting_horizon]
        )

        poisson_dev_out <- poisson_dev_out -2*(models_ll-saturated_model_ll)

      }

      poisson_dev_out_g <- -2*(models_ll-saturated_model_ll)

      l[[paste0('oos_deviance_',i)]] <- poisson_dev_out_g
      l[[paste0('error_',i)]] <- relative_mse_g

      }

    l[['oos_deviance']] <- poisson_dev_out
    l[['error']] <- relative_mse / (denominator)


    out[[model_ix]] <- l

  }


  return(out)



}


model_assessment_superpop <- function(model_fit_and_prediction,
                             N_groups,
                             years_fit,
                             forecasting_horizon,
                             bias=15,
                             subset_rates=NULL) {
  out <- list()

  relative_mse <- 0
  denominator <- 0
  poisson_dev_out <- 0
  tmp_actual_data <- model_fit_and_prediction[['actual_data']]


  for (model_ix in c("lc", "apc", "rh")) {
    tmp_model <- model_fit_and_prediction[[model_ix]]

    l <- list()

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i = "superpop"

    tmp_group_string_actual <- paste0("muxt_actual_", i)
    tmp_group_string_fitted <- paste0("muhat", i)

    if(!is.null(subset_rates)){

      subsetting_condition <- rep(FALSE, length(tmp_actual_data[[tmp_group_string_actual]]))

      subsetting_condition[subset_rates-bias+1] <- TRUE

    }else{

      subsetting_condition <- rep(TRUE, length(tmp_actual_data[[tmp_group_string_actual]]))

    }

    condition_selection <- !is.na(tmp_actual_data[[tmp_group_string_actual]]) &
      tmp_actual_data[[tmp_group_string_actual]] > 0 &subsetting_condition

    denominator <- denominator + sum(condition_selection)

    relative_mse_g <- (sum(((abs(
      tmp_model[[tmp_group_string_fitted]] - tmp_actual_data[[tmp_group_string_actual]]
    )) / tmp_actual_data[[tmp_group_string_actual]])[condition_selection]))/sum(condition_selection)


    models_ll <- -poisson_nll(
      occurrence = tmp_actual_data$full_pp_data$datahat_tot$Dxt[condition_selection, length(years_fit) +
                                                                  forecasting_horizon],
      exposure = tmp_actual_data$full_pp_data$datahat_tot$Ext[condition_selection, length(years_fit) +
                                                                forecasting_horizon],
      muxt_hat = tmp_model[[tmp_group_string_fitted]][condition_selection]
    )

    saturated_model_ll <- -poisson_nll(
      occurrence = tmp_actual_data$full_pp_data$datahat_tot$Dxt[condition_selection, length(years_fit) +
                                                                  forecasting_horizon],
      exposure = tmp_actual_data$full_pp_data$datahat_tot$Ext[condition_selection, length(years_fit) +
                                                                forecasting_horizon],
      muxt_hat = tmp_actual_data$full_pp_data$datahat_tot$Dxt[condition_selection, length(years_fit) +
                                                                forecasting_horizon]/tmp_actual_data$full_pp_data$datahat_tot$Ext[condition_selection, length(years_fit) +
                                                                                                                                    forecasting_horizon]
    )

    poisson_dev_out_superpop <-  -2*(models_ll-saturated_model_ll)


    l[[paste0('oos_deviance')]] <- poisson_dev_out_superpop
    l[[paste0('error')]] <- relative_mse_g
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for(i in 1:N_groups){

      l[[paste0('oos_deviance_',i)]] <- NA
      l[[paste0('error_',i)]] <- NA
    }

    out[[model_ix]] <- l

  }


  return(out)



}





model_assessment_full_mle <- function(model_fit_and_prediction,
                             N_groups,
                             years_fit,
                             forecasting_horizon,
                             bias=15,
                             subset_rates=NULL) {
  out <- list()

  relative_mse <- 0
  denominator <- 0
  poisson_dev_out <- 0
  tmp_actual_data <- model_fit_and_prediction[['actual_data']]


  for (model_ix in c("lc", "apc", "rh")) {
    tmp_model <- model_fit_and_prediction[[model_ix]]

    l <- list()

    for (i in 1:N_groups) {
      tmp_group_string_actual <- paste0("muxt_actual_", i)
      tmp_group_string_fitted <- paste0("muhat", i, "_MLE")

      if(!is.null(subset_rates)){

        subsetting_condition <- rep(FALSE, length(tmp_actual_data[[tmp_group_string_actual]]))

        subsetting_condition[subset_rates-bias+1] <- TRUE

      }else{

        subsetting_condition <- rep(TRUE, length(tmp_actual_data[[tmp_group_string_actual]]))

      }

      # browser()
      condition_selection <- !is.na(tmp_actual_data[[tmp_group_string_actual]]) &
        tmp_actual_data[[tmp_group_string_actual]] > 0 &subsetting_condition


      denominator <- denominator + sum(condition_selection)


      relative_mse <- relative_mse + sum(((abs(
        tmp_model[[tmp_group_string_fitted]] - tmp_actual_data[[tmp_group_string_actual]]
      )) / tmp_actual_data[[tmp_group_string_actual]])[condition_selection])

      relative_mse_g <- (sum(((abs(
        tmp_model[[tmp_group_string_fitted]] - tmp_actual_data[[tmp_group_string_actual]]
      )) / tmp_actual_data[[tmp_group_string_actual]])[condition_selection]))/sum(condition_selection)


      if (i < 2) {
        models_ll <- -poisson_nll(
          occurrence = tmp_actual_data$full_pp_data$datahat$Dxt[condition_selection, length(years_fit) +
                                                                  forecasting_horizon],
          exposure = tmp_actual_data$full_pp_data$datahat$Ext[condition_selection, length(years_fit) +
                                                                forecasting_horizon],
          muxt_hat = tmp_model[[tmp_group_string_fitted]][condition_selection]
        )

        saturated_model_ll <- -poisson_nll(
          occurrence = tmp_actual_data$full_pp_data$datahat$Dxt[condition_selection, length(years_fit) +
                                                                  forecasting_horizon],
          exposure = tmp_actual_data$full_pp_data$datahat$Ext[condition_selection, length(years_fit) +
                                                                forecasting_horizon],
          muxt_hat = tmp_actual_data$full_pp_data$datahat$Dxt[condition_selection, length(years_fit) +
                                                                forecasting_horizon]/tmp_actual_data$full_pp_data$datahat$Ext[condition_selection, length(years_fit) +
                                                                                                                                forecasting_horizon]
        )

        poisson_dev_out <- poisson_dev_out  -2*(models_ll-saturated_model_ll)

      } else{
        models_ll <- -poisson_nll(
          occurrence = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i - 1]][['Dxt']][condition_selection, length(years_fit) +
                                                                                                forecasting_horizon] ,
          exposure = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i -
                                                                             1]][['Ext']][condition_selection, length(years_fit) +
                                                                                            forecasting_horizon],
          muxt_hat = tmp_model[[tmp_group_string_fitted]][condition_selection]
        )

        saturated_model_ll <- -poisson_nll(
          occurrence = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i - 1]][['Dxt']][condition_selection, length(years_fit) +
                                                                                                forecasting_horizon] ,
          exposure = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i -
                                                                             1]][['Ext']][condition_selection, length(years_fit) +
                                                                                            forecasting_horizon],
          muxt_hat = tmp_actual_data$full_pp_data$list_of_extra_exposures[[i - 1]][['Dxt']][condition_selection, length(years_fit) +
                                                                                              forecasting_horizon]/tmp_actual_data$full_pp_data$list_of_extra_exposures[[i -
                                                                                                                                                                           1]][['Ext']][condition_selection, length(years_fit) +
                                                                                                                                                                                          forecasting_horizon]
        )

        poisson_dev_out <- poisson_dev_out -2*(models_ll-saturated_model_ll)

      }


      poisson_dev_out_g <- -2*(models_ll-saturated_model_ll)

      l[[paste0('oos_deviance_',i)]] <- poisson_dev_out_g
      l[[paste0('error_',i)]] <- relative_mse_g

    }

    l[['oos_deviance']] <- poisson_dev_out
    l[['error']] <- relative_mse / (denominator)

    out[[model_ix]] <- l

  }


  return(out)



}

# Plots ----

weights_plotter <- function(credibility_model,
                            N_groups,
                            ages_fit,
                            predictor = "lc",
                            ages_breaks=c(50, 60, 70, 80, 90)) {
  tmp <- credibility_model[[predictor]]


  zeds <- c()
  names_to_add <- c()


  for (i in 1:N_groups) {
    weights_to_subset <- paste0("Z_", i)
    names_to_add <- c(names_to_add, weights_to_subset)

    zeds <- c(zeds, 1 - tmp[[weights_to_subset]])

  }

  dt_Zs <- data.frame(
    Zs = zeds,
    label = factor(rep(names_to_add, each = length(ages_fit)),
                   levels = c("Z_2","Z_3","Z_1")),
    ages.code = rep(ages_fit, N_groups)
  )


  text_size = 28


  dt_Zs %>%
    ggplot(aes(x = ages.code, y = Zs)) +
    geom_point(aes(colour = label), size = 3, alpha = .7) +
    theme_bw() +
    scale_x_continuous(breaks = ages_breaks) +

    theme(
      text = element_text(size = text_size),
      legend.position = "inside",                          # must now be "inside"
      # legend.position.inside = c(0.85, 0.85),              # numeric position
      legend.justification = c("right", "bottom"),
      legend.background = element_blank(),   # removes legend box
      legend.key = element_blank()
    ) +
    ylab("") +
    xlab("") +
    labs(color = "")




}


all_thetas_plotter <- function(credibility_model,
                               N_groups,
                               ages_fit,
                               predictor = "lc") {
  tmp <- credibility_model[[predictor]]


  cis <- c()
  names_to_add <- c()


  for (i in 1:N_groups) {
    cs_to_subset <- paste0("C", i)
    names_to_add <- c(names_to_add, cs_to_subset)

    cis <- c(cis, tmp[[cs_to_subset]])

  }

  dt_cis <- data.frame(
    cis = cis,
    label = factor(rep(names_to_add, each = length(ages_fit)),
                   levels=c('C2', 'C3', 'C1')),
    ages.code = rep(ages_fit, N_groups)
  )



  text_size = 28


  dt_cis %>%
    ggplot(aes(x = ages.code, y = cis)) +
    geom_point(aes(colour = label), size = 3, alpha = .7) +
    ylim(0.5, 2.1) +
    theme_bw() +
    scale_x_continuous(breaks = (1:9)*10) +
    theme(
      text = element_text(size = 28),
      legend.position = "inside",                          # must now be "inside"
      # legend.position.inside = c(0.85, 0.85),              # numeric position
      legend.justification = c("right", "top"),
      legend.key.height = unit(1.5, "lines"),
      legend.background = element_blank(),   # removes legend box
      legend.key = element_blank()
    ) +
    ylab("") +
    xlab("") +
    labs(color = "")




}



exposure_plotter <- function(data_pp,
                           predictor = "lc") {


  E1 <- apply(data_pp$datahat$Ext, 1, sum, na.rm = T)
  E2 <-  apply(data_pp$list_of_extra_exposures[[1]]$Ext, 1, sum, na.rm =
                 T)
  E3 <-  apply(data_pp$list_of_extra_exposures[[2]]$Ext, 1, sum, na.rm =
                 T)


  dt <- data.frame(
    exposure = c(E1, E2, E3),
    Age = rep(names(E1), 3),
    group = rep(c("p", "p1", "p2"), each = length(names(E1)))
  )

  p <- dt %>%
    mutate(
      group=recode(group,
                   p = "Group 0",
                   p1 = "Group 1",
                   p2 = "Group 2")) %>%
    ggplot(aes(x=as.integer(Age), y=exposure)) +
    geom_line(aes(colour = group), size = 1.5)+
    theme_bw()+
    theme(text = element_text(size = 28),
          legend.position = "top") +
    scale_color_manual(values = c("Group 0" = "#a71429", "Group 1" = "#4169E1", "Group 2" = "#2E8B57"))+
    ylab("") +
    xlab("") +
    labs(color = "")

  return(p)


}



variance_plotter <- function(credibility_model,
                             predictor = "lc",
                             ages_fit) {

  # browser()
  v1 <- credibility_model[[predictor]]$varthetax_1
  v2 <-  credibility_model[[predictor]]$varthetax_2
  v3 <-  credibility_model[[predictor]]$varthetax_3


  dt <- data.frame(
    variance = c(v1,v2,v3),
    Age = rep(as.character(ages_fit), 3),
    group = rep(c("p", "p1", "p2"), each = length(as.character(ages_fit))))


  p <- dt %>%
    mutate(
      group = recode(group,

                     p1 = "Group 1",
                     p2 = "Group 2",
                     p = "Group 3")
    ) %>%
    ggplot(aes(x = as.integer(Age), y = (variance),
               color = interaction(group), shape = interaction(group))) +
    geom_point(size = 3) +
    theme_bw() +
    theme(
      text = element_text(size = 28),
      legend.position = "inside",                          # must now be "inside"
      # legend.position.inside = c(0.85, 0.85),              # numeric position
      legend.justification = c("right", "top"),
      legend.background = element_blank(),   # removes legend box
      legend.key = element_blank()
    )+
    scale_color_manual(values = c("Group 1" = "#4169E1", "Group 2" = "#2E8B57","Group 3" = "#a71429"),
                       labels = c(expression(hat(Var)(Theta[x]^1)), expression(hat(Var)(Theta[x]^2)),expression(hat(Var)(Theta[x]^3)))
    ) +
    scale_shape_manual(values = c("Group 1" = 16, "Group 2" = 17, "Group 3" = 18),
                       labels = c(expression(hat(Var)(Theta[x]^1)), expression(hat(Var)(Theta[x]^2)),expression(hat(Var)(Theta[x]^3)))) +
    labs(color = "Group", shape = "Group") +  # Customize the legend titles
    ylab("") +
    xlab("") +
    # ylim(-0.0001, 1)+
    labs(color = "",shape="")
    guides(color = guide_legend(override.aes = list(shape = c(16, 17, 18))))  # Ensures the legend reflects both shape and color


  return(p)


}


thetas_plotter <- function(credibility_model,
                           subgroup,
                           ages_fit,
                           years_fit,
                           mortality_models_fit,
                           predictor = "lc",
                           chosen_age = 70) {


  tmp <- credibility_model[[predictor]]

  tmp_fit <- mortality_models_fit[[predictor]]


  mortality_model_forecast <- forecast(tmp_fit, h = 1 #doesnt really matter here. I want the fitted values
                                       )

 muxt_hat <- mortality_model_forecast$fitted


 chosen_age = as.character(chosen_age)

 selection_ratios <- paste0("Fxt_", subgroup)



 Cs_frame <- data.frame(y = (tmp[[selection_ratios]] * 1 / muxt_hat)[chosen_age, ], x =
                          years_fit) %>%
   filter(y > 0)
 selection_ratios <- paste0("C", subgroup)
 C1 <- tmp[[selection_ratios]]
 selection_ratios <- paste0("varthetax_", subgroup)
 varthetax_1 <- tmp[[selection_ratios]]

 names(varthetax_1) <- names(C1)

 text_size = 28

 Cs_frame %>%
   ggplot(aes(x = x, y = y)) +
   geom_point(size = 3) +
   geom_hline(
     yintercept = C1[chosen_age],
     linetype = "dotted",
     linewidth = 1.2,
     color = "#a71429"
   ) +
   geom_hline(
     yintercept = C1[chosen_age] + sqrt(varthetax_1[chosen_age]),
     linewidth = 1.2,
     color = "#4169E1",
     linetype = "dotted"
   ) +
   geom_hline(
     yintercept = C1[chosen_age] - sqrt(varthetax_1[chosen_age]),
     linewidth = 1.2,
     color = "#4169E1",
     linetype = "dotted"
   ) +
   theme_bw() +
   scale_x_continuous(breaks = c(10, 20, 30, 40, 50)) +
   theme(text = element_text(size = text_size),
         legend.position = "top") +
   ylab("") +
   xlab("") +
   labs(color = "")

}




binomial_simulator <- function(starting_exposure, probabilities) {
  starting_period_lives <- starting_exposure
  starting_period_events <- NULL

  for (ix in 1:length(probabilities)) {
    events <- rbinom(1, size = starting_period_lives[ix], p = probabilities[ix])

    starting_period_lives <- c(starting_period_lives, starting_period_lives[ix] -
                                 events)

    starting_period_events <- c(starting_period_events, events)

  }


  l <- list(Lx = starting_period_lives[-length(starting_period_lives)], Dx =
              starting_period_events)

  return(l)


}


average_exposure_and_deaths <- function(raw_exposure) {
  out <- raw_exposure
  events <- NULL



  # One would need some adjustments for older ages and age 0. Since we
  # do not care about those we simply disgregard them.
  for (ix in 1:(length(raw_exposure)) - 1) {


    out[ix] <- raw_exposure[ix] / 3 + raw_exposure[ix + 1] / 6 + raw_exposure[ix] /
      6 + raw_exposure[ix + 1] / 3

    events[ix] <- (raw_exposure[ix] - raw_exposure[ix + 1]) / 2
  }



  events <- c(events,last(raw_exposure))
    l<- list(Ext=as.numeric(out),
             Dxt=as.numeric(events))
  return(l)

}


plot_predicted_event_trends <- function(mortality_models_fit,
                                        data_pp,
                                        # subgroup=1,
                                        model_option ="apc"){
  # browser()
  mortality_model_fit <- mortality_models_fit[[model_option]]


  mortality_model_forecast <- forecast(mortality_model_fit,
                                       h = 1)


  # if (subgroup != 1) {
  #   tmp_exp <- data_pp$list_of_extra_exposures[[subgroup-1]]$Ext
  #
  # } else{
  #   tmp_exp <- data_pp$datahat$Ext
  # }

  tmp_exp <- data_pp$datahat$Ext

  tmp1 <- data.frame(y=apply(tmp_exp*mortality_model_forecast$fitted,1, sum,na.rm=T),
                   x=rownames(tmp_exp),
                   group="Group 3")

  tmp_exp <- data_pp$list_of_extra_exposures[[1]]$Ext

  tmp2 <- data.frame(y=apply(tmp_exp*mortality_model_forecast$fitted,1, sum,na.rm=T),
                     x=rownames(tmp_exp),
                     group="Group 1")

  tmp_exp <- data_pp$list_of_extra_exposures[[2]]$Ext

  tmp3 <- data.frame(y=apply(tmp_exp*mortality_model_forecast$fitted,1, sum,na.rm=T),
                     x=rownames(tmp_exp),
                     group="Group 2")

  dt <- rbind(tmp1,tmp2,tmp3)

  p <- dt %>%
    ggplot(aes(x = as.numeric(x), y = log(y), color=group)) +
    geom_line(size = 2, alpha = .7) +
    theme_bw() +
    # scale_x_continuous(breaks = ages_breaks) +
    scale_color_manual(values = c("Group 1" = "#4169E1", "Group 2" = "#2E8B57","Group 3" = "#a71429"),
                       labels = c(expression(D[x,.]^1), expression(D[x,.]^2), expression(D[x, .]^3)))+
    theme(
      text = element_text(size = 28),
      legend.position = "inside",                          # must now be "inside"
      # legend.position.inside = c(0.85, 0.85),              # numeric position
      legend.justification = c("right", "bottom"),
      legend.background = element_blank(),   # removes legend box
      legend.key = element_blank()
    ) +
    scale_fill_discrete(guide="none")+
    ylab("") +
    # ylab(expression(log(hat(D)[x* "."]^i ))) +
    xlab("") +
    labs(color = "")

  return(p)

}

fill_missing_rows <- function(df) {
  # First, create a vector of ages from 0 to 110
  all_ages <- 0:110

  # Initialize an empty list to store results
  df_filled <- list()

  # For each cohort, ensure all ages are represented in the data
  for (cohort in unique(df$Cohort)) {
    # Filter data for the current cohort
    cohort_data <- df %>% filter(Cohort == cohort)

    # Get the set of ages already present for the current cohort
    present_ages <- unique(cohort_data$Age)

    # Find missing ages by subtracting present ages from all ages (0 to 110)
    missing_ages <- setdiff(all_ages, present_ages)

    # If there are missing ages, create new rows for the missing ages
    if (length(missing_ages) > 0) {
      # Create new rows for the missing ages (using data from Cohort 0 for 'qx' values)
      missing_rows <- data.frame(
        Age = missing_ages,
        Period = missing_ages + cohort,  # Correct Period as Cohort + Age
        Cohort = cohort,
        qx = df$qx[df$Cohort == 0 & df$Age %in% missing_ages]  # Match missing age rows from Cohort 0
      )

      # Combine the existing cohort rows with the missing rows (do not duplicate)
      df_filled[[length(df_filled) + 1]] <- bind_rows(
        cohort_data,  # Keep all existing rows in cohort_data
        missing_rows  # Add only the missing rows
      )
    } else {
      # If no ages are missing, just add the current cohort's data to the list
      df_filled[[length(df_filled) + 1]] <- cohort_data
    }
  }

  # Combine all the dataframes in the list into one dataframe
  df_filled <- bind_rows(df_filled)

  return(df_filled)
}



k_years_smoothing <- function(x,k=5){

  ll <- 1:length(x)

  out <- sapply(ll,average_smoother,x=x,k=k)

  return(out)

}


average_smoother <- function(pos,x,k){

  out <- x[max(pos-((k-1)/2),1):min(pos+((k-1)/2),length(x))]
  return(mean(out))
}


plot_smoothed_unsmoothed  <- function(data,
                                                                           data_pp,
                                                                           N_groups,
                                                                           mortality_models_fit,
                                                                           years_fit,
                                                                           ages_fit,
                                                                           bias=70,
                                                                           scenario = 2,
                                      model_2_plot="lc",
                                                                           smoothing_years=5,
                                                                           forecasting_horizon) {



  out <- list()
  observed_rates <- list()

  data_pp_2 <- data_preprocessing_scenario_2(
    data = data,
    N_groups = N_groups,
    bias = bias,
    ages_fit = ages_fit,
    years_fit = c(years_fit, (last(years_fit) + 1):(last(years_fit) + forecasting_horizon))
    # scenario = scenario
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

    # browser()

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

    tmp_occ <- data_pp$datahat$Dxt[as.character(ages_fit), as.character(years_fit)]
    tmp_exp <- data_pp$datahat$Ext[as.character(ages_fit), as.character(years_fit)]

    C1 <- apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp * muxt_hat, 1, sum, na.rm =
                                                      T)
    # browser()
    l[['C1_unsmoothed']] <- C1
    fit1 <- rpart(C1~ages_fit)
    C1 <- predict(fit1, data.frame(ages_fit))


    # C1 <- k_years_smoothing(C1, k= smoothing_years)

    Fxt_1 <- tmp_occ / tmp_exp

    tmpvar1 <- ((Fxt_1 - muxt_hat)^2 - (muxt_hat / tmp_exp))
    # tmpvar1[tmpvar1 == 1] <- NA
    varthetax_1 <- pmax(0,
                        apply(tmpvar1, 1, sum, na.rm = T) / apply(muxt_hat^2, 1, sum, na.rm = T),
                        na.rm = T)

    # varthetax_1 <- k_years_smoothing(varthetax_1, k= smoothing_years)
    l[['varthetax_1_unsmoothed']] <- varthetax_1
    fit1 <- rpart(varthetax_1~ages_fit)
    varthetax_1 <- predict(fit1, data.frame(ages_fit))
    l[['varthetax_1_smoothed']] <- varthetax_1

    Z_1 <- 1 / (1 + varthetax_1 * apply(tmp_exp * muxt_hat, 1, sum, na.rm =
                                          T))

    subgroups_data <- data_pp$list_of_extra_exposures


    l[['C1_smoothed']] <- C1



    # term1 <- 2*(muxt_hat_predicted^2)*(varthetax_1)+muxt_hat_predicted^2
    #
    # term2 <- varthetax_1*(apply((tmp_exp^2) *( muxt_hat)^2, 1, sum, na.rm =
    #                   T)/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
    #                     T)^2)+(1/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
    #                                     T))
    #
    # l[['msep_1']] <- term1+term2


    # observed_rates[['muxt_actual_1']] <- data_pp_2$datahat$Dxt[, (length(years_fit) +
    #                                                                 1):(length(years_fit) + forecasting_horizon)] / data_pp_2$datahat$Ext[, (length(years_fit) +
    #                                                                                                                                            1):(length(years_fit) + forecasting_horizon)]
    observed_rates[['muxt_actual_1']] <- data_pp_2$datahat$Dxt[, (length(years_fit) + forecasting_horizon)] / data_pp_2$datahat$Ext[, (length(years_fit) + forecasting_horizon)]

    for (i in 2:N_groups) {
      tmp_occ <- subgroups_data[[i - 1]]$Dxt[as.character(ages_fit), as.character(years_fit)]
      tmp_exp <- subgroups_data[[i - 1]]$Ext[as.character(ages_fit), as.character(years_fit)]

      assign(
        paste0("C", i),
        apply(tmp_occ, 1, sum, na.rm = T) / apply(tmp_exp *
                                                    muxt_hat, 1, sum, na.rm = T)
      )


      l[[paste0('C',i,'_unsmoothed')]] <- get(paste0("C", i))
      fit1 <- rpart(get(paste0("C", i))~ages_fit)
      assign(paste0("C", i), predict(fit1, data.frame(ages_fit)))

      assign(paste0("Fxt_", i), tmp_occ / tmp_exp)

      assign(paste0("tmpvar", i), ((get(
        paste0("Fxt_", i)
      ) - muxt_hat)^2 - (muxt_hat / tmp_exp)))

      tmp <- get(paste0("tmpvar", i))
      # tmp[tmp == 1] = NA
      assign(paste0("tmpvar", i), tmp)
      assign(paste0("varthetax_", i), pmax(
        0,
        apply(get(paste0("tmpvar", i)), 1, mean, na.rm =
                T) / apply(muxt_hat^2, 1, mean, na.rm =
                             T),
        na.rm = T
      ))

      l[[(paste0("varthetax_",i,"_unsmoothed"))]] <- get(paste0("varthetax_", i))
      fit1 <- rpart(get(paste0("varthetax_", i))~ages_fit)
      assign(paste0("varthetax_", i), predict(fit1, data.frame(ages_fit)))
      l[[(paste0("varthetax_",i,"_smoothed"))]] <- get(paste0("varthetax_", i))


      l[[paste0('C',i,'_smoothed')]] <- get(paste0("C", i))

      assign(paste0("Z_", i), 1 / (1 + get(paste0(
        "varthetax_", i
      )) *
        apply(tmp_exp * muxt_hat, 1, sum, na.rm = T)))




      # l[[paste0("C", i)]] <- get(paste0("C", i))
      # l[[paste0("Fxt_", i)]] <- get(paste0("Fxt_", i))
      # l[[paste0("Z_", i)]] <- get(paste0("Z_", i))
      # l[[paste0("varthetax_", i)]] <- get(paste0("varthetax_", i))

      # l[[paste0("muhat", i)]] <- get(paste0("Z_", i)) * muxt_hat_predicted + (1 - get(paste0("Z_", i))) * get(paste0("C", i)) * muxt_hat_predicted
      # l[[paste0("muhat", i,"_MLE")]] <- get(paste0("C", i)) * muxt_hat_predicted


      # term1 <- 2*(muxt_hat_predicted^2)*(get(paste0("varthetax_", i)))+muxt_hat_predicted^2
      #
      # term2 <- get(paste0("varthetax_", i))*(apply((tmp_exp^2) *( muxt_hat)^2, 1, sum, na.rm =
      #                               T)/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
      #                                           T)^2)+(1/ apply(tmp_exp * muxt_hat, 1, sum, na.rm =
      #                                                             T))
      #
      # l[[paste0('msep_',i)]] <- term1+term2

      subgroups_forecast <- data_pp_2$list_of_extra_exposures[[i - 1]]

      # observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) +
      #                                                                            1):(length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) +
      #                                                                                                                                                        1):(length(years_fit) + forecasting_horizon)]
      #
      observed_rates[[paste0('muxt_actual_', i)]] <- subgroups_forecast$Dxt[, (length(years_fit) + forecasting_horizon)] / subgroups_forecast$Ext[, (length(years_fit) + forecasting_horizon)]

    }

    # observed_rates[['full_pp_data']] <- data_pp_2
    # out[['actual_data']] <- observed_rates
    out[[model_option]] <- l

  }



  selected_out <- out[[model_2_plot]]



  dt_2_plot <- data.frame(

    population = rep(c(1,2,3), each=length(ages_fit)),
    theta_unsmoothed= c(selected_out$C1_unsmoothed,
                     selected_out$C2_unsmoothed,
                     selected_out$C3_unsmoothed),
    theta_smoothed= c(selected_out$C1_smoothed,
                     selected_out$C2_smoothed,
                     selected_out$C3_smoothed),
    variance_unsmoothed= c(selected_out$varthetax_1_unsmoothed,
                     selected_out$varthetax_2_unsmoothed,
                     selected_out$varthetax_3_unsmoothed),
    variance_smoothed= c(selected_out$varthetax_1_smoothed,
                         selected_out$varthetax_2_smoothed,
                         selected_out$varthetax_3_smoothed)

  ) %>%
    mutate(    group_variance=cumsum(c(TRUE, diff(variance_smoothed) != 0)),
               group_theta=cumsum(c(TRUE, diff(theta_smoothed) != 0)),
               age_code= rep(as.character(ages_fit),3))

  dt_1 <- dt_2_plot %>%
    select(theta_unsmoothed,
           age_code,
           population)
  dt_2 <- dt_2_plot %>%
    select(theta_smoothed,
           group_theta,
           age_code,
           population)


  p_theta <- ggplot(dt_1, aes(x = age_code, y = theta_unsmoothed,color=as.character(population))) +
    geom_point(size = 2, alpha = .3) +
    geom_line(data=dt_2,
              aes(x = age_code, y = theta_smoothed, group = group_theta,color=as.character(population)),
              size=2,
              alpha=.7) +
    ylim(0.2, 2.3) +
    theme_bw() +
    scale_x_discrete(breaks = c(20,40,60,80,100)) +
    theme(text = element_text(size = 28),
          legend.position = "top") +
    ylab("") +
    xlab("") +
    labs(color = "")+
    scale_color_manual(
      values = c("1" = "#a71429","2" = "#4169E1", "3" = "#006400"),
      labels = c(expression(hat(theta)[x]^0), expression(hat(theta)[x]^1), expression(hat(theta)[x]^2))
    )


  dt_1 <- dt_2_plot %>%
    select(variance_unsmoothed,
           age_code,
           population)
  dt_2 <- dt_2_plot %>%
    select(variance_smoothed,
           group_variance,
           age_code,
           population)

  p_var <- ggplot(dt_1, aes(x = age_code, y = variance_unsmoothed,color=as.character(population))) +
    geom_point(size = 2, alpha = .3) +
    geom_line(data=dt_2,
              aes(x = age_code, y = variance_smoothed, group = group_variance,color=as.character(population)),
              size=2,
              alpha=.7) +
    ylim(-.0000001, 2.3) +
    theme_bw() +
    scale_x_discrete(breaks = c(20,40,60,80,100)) +
    theme(text = element_text(size = 28),
          legend.position = "top") +
    # ylim(0,.7)+
    ylab("") +
    xlab("") +
    labs(color = "")+
    scale_color_manual(
      values = c("1" = "#a71429","2" = "#4169E1", "3" = "#006400"),
      labels = c(expression(hat(Var)(Theta[x]^0)), expression(hat(Var)(Theta[x]^1)), expression(hat(Var)(Theta[x]^2)))
    )

  return(list(variance_plot=p_var,
              theta_plot= p_theta))
}

