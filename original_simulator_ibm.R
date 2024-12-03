rm(list = ls())

library(IBMPopSim)
library(demography)

{
  set.seed(71)
  N <- 1000
  pop_df <- data.frame(
    "birth" = sample(-65, size = N, replace = T),
    "death" = rep(NA, N),
    "risk_cls" = sample(
      c(1, 2),
      size =
        N,
      prob =
        c(0.5, 0.5),
      replace =
        T
    )
  ) # rep(1:2,each= N/2)
  pop_init <- population(pop_df, entry = TRUE, out = TRUE)
  LC <- lc()
  ages.fit <- 65:100
  years.fit <- 1950:2018
  
  ICE.mx <- hmd.mx(country = "SWE",
                   "gabriele.pittarello@uniroma1.it",
                   "pm869fqW8nxozS8!")
  D5data <- EWStMoMoMale <- StMoMoData(ICE.mx, series = 'male', type = "central")
  
  LCfitMale <- fit(LC,
                   data = D5data,
                   ages.fit = ages.fit,
                   years.fit = years.fit)
  
  # EWStMoMoMale <- StMoMoData(EWdata_hmd, series = "male")
  #Fitting
  LC <- lc()
  ages.fit <- 65:100
  years.fit <- 1950:2018
  LCfitMale <- fit(LC,
                   data = EWStMoMoMale,
                   ages.fit = ages.fit,
                   years.fit = years.fit)
  ## StMoMo: Start fitting with gnm
  ## Initialising
  ## Running start-up iterations..
  ## Running main iterations.....
  ## Done
  ## StMoMo: Finish fitting with gnm
  t <- 30
  LCforecastMale <- forecast(LCfitMale, h = t)
  
  d_k <- apply(LCforecastMale$rates, 2, function(x)
    stepfun(66:100, x))
  breaks <- 1:29
  death_male <- piecewise_xy(breaks, d_k)
  death_male(10, 65) # Death rate at time t=10 (years 2027) and age 65.
  params <- list("death_male" = death_male, "alpha" = c(1.3, 0.8))
  
  
  
  params$mu <- c(0.001, 0.06)
  params$lambda <- 900 # Entry events
  
  
  death_event <- mk_event_individual(type = "death", intensity_code = "result = alpha[I.risk_cls-1] * death_male(t,age(I, t));")
  
  
  exit_event <- mk_event_individual(type = "exit", intensity = "result = mu[I.risk_cls-1]; ")
  
  
  params$p <- 0.5
  entry_event <- mk_event_poisson(
    type = "entry",
    intensity = "lambda",
    kernel_code = "if (CUnif()<p)
newI.risk_cls =1;
else
newI.risk_cls= 2;
double a = CUnif(65,70);
newI.set_age(a,t);
newI.out = false;"
  )
  
  
  
  model <- mk_model(
    characteristics = get_characteristics(pop_init),
    # Characteristics names and types
    events = list(death_event, entry_event, exit_event),
    # Events list
    parameters = params # Model parameters
  )
  
  
  death_max <- max(sapply(d_k, function(x) {
    max(x)
  }))
  
  
  sim_out <- popsim(
    model = model,
    initial_population = pop_init,
    events_bounds = c(
      'death' = death_max,
      "entry" = params$lambda,
      "exit" = max(params$mu)
    ),
    parameters = params,
    time = t,
    age_max = 110,
    multithreading = TRUE
  )
}


ages.fit <- 65:80
years.fit <- 0:65

N_groups <- 3

# find deaths and exposures data

D1 <- death_table(
  sim_out$population[sim_out$population$risk_cls == 1, ],
  ages = c(ages.fit, last(ages.fit)),
  period = c(years.fit, last(years.fit))
)

D2 <- death_table(
  sim_out$population[sim_out$population$risk_cls == 2, ],
  ages = c(ages.fit, last(ages.fit)),
  period = c(years.fit, last(years.fit))
)

E1 <- exposure_table(
  sim_out$population[sim_out$population$risk_cls == 1, ],
  ages = c(ages.fit, last(ages.fit)),
  period = c(years.fit, last(years.fit))
)

E2 <- exposure_table(
  sim_out$population[sim_out$population$risk_cls == 2, ],
  ages = c(ages.fit, last(ages.fit)),
  period = c(years.fit, last(years.fit))
)


# create a super-population

Dtot <- EWStMoMoMale$Dxt[as.character(ages.fit), as.character(1950:(1950 +
                                                                      max(years.fit)))]
Etot <- EWStMoMoMale$Ext[as.character(ages.fit), as.character(1950:(1950 +
                                                                      max(years.fit)))]

colnames(Etot) <- colnames(Dtot) <- years.fit[1:(length(years.fit))]


datahat <- structure(
  list(
    Dxt = Dtot - D1 - D2 ,
    Ext = Etot - E1 - E2,
    ages = ages.fit,
    years = years.fit,
    type = 'central',
    series = 'female',
    label = 'total'
  ),
  class = "StMoMoData"
)


mortality_model_lc <- lc(link = "log")
model_option <- "lc"
assign("mortality_model", get(paste0("mortality_model_", model_option)))


list_of_extra_exposures <- list()

for (i in 2:N_groups) {
  # i=3
  tmp_list <- list()
  
  tmp_list[['Dxt']] <- get(paste0("D", i - 1))
  tmp_list[['Ext']] <- get(paste0("E", i - 1))
  
  list_of_extra_exposures[[i - 1]] <- tmp_list
  
}


mortality_model_fit <- fit(
  mortality_model,
  data = datahat,
  years.fit = years.fit,
  ages.fit = ages.fit,
  list_of_extra_exposures = list_of_extra_exposures
)

forecasting_horizon <- 1

year.predict <- max(years.fit) + forecasting_horizon

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



C1 <- apply(D1, 1, sum) / apply(E1 * muxt_hat, 1, sum)

C2 <- apply(D2, 1, sum) / apply(E2 * muxt_hat, 1, sum)

Fxt_1 <- D1 / E1
Fxt_2 <- D2 / E2


tmpvar1 <- ((Fxt_1 - muxt_hat) ^ 2) / (muxt_hat ^ 2)
tmpvar2 <- ((Fxt_2 - muxt_hat) ^ 2) / (muxt_hat ^ 2)
tmpvar1[tmpvar1 == 1 | tmpvar1 == 0] <- NA
tmpvar2[tmpvar2 == 1 | tmpvar2 == 0] <- NA

varthetax_1 <- apply(tmpvar1, 1, mean, na.rm = T)
varthetax_2 <- apply(tmpvar2, 1, mean, na.rm = T)


Z_1 <- 1 / (1 + varthetax_1 * apply(E1 * muxt_hat, 1, sum))

Z_2 <- 1 / (1 + varthetax_2 * apply(E2 * muxt_hat, 1, sum))

1 - Z_1
1 - Z_2
varthetax_1
varthetax_2
