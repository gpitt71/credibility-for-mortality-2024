# Libraries  ----
rm(list=ls())
Sys.setenv(LANGUAGE = "en")
library(IBMPopSim)
library(dplyr)
library(StMoMo)
library(ggplot2)

# Group and age effect -----


data_generator_0 <- function(age_sample_size = 100000,
                             super_pop_share = 0.8,
                             seed = 1996){
  
  N <- age_sample_size
  
  
  Q1 <- super_pop_share
  N1 <- round(N*Q1)
  ages_1 <- -80:0
  tmp_df <- data.frame("birth"=rep(ages_1,each=N),
                       "death"=rep(rep(NA,N),each=length(ages_1)),
                       "risk_cls"= 1)
  
  
  ages <- -65:-50
  age_classes <- length(ages)
  
  pop_df <- data.frame("birth"=rep(ages,each=N-N1),
                       "death"=rep(rep(NA,N-N1),each=age_classes),
                       "risk_cls"= sample(c(2,3), 
                                          size=(N-N1)*age_classes,
                                          replace=T,
                                          prob=c(.8,.2)))
  

  
  superpop_df <- population(tmp_df, entry=TRUE, out=TRUE)
  
  subpop_df <- population(pop_df, entry=TRUE, out=TRUE)
  
  EWStMoMoMale <- StMoMoData(EWdata_hmd, series = "male")
  
  
  #Fitting
  LC <- lc()
  ages.fit <- 50:100
  years.fit <- 1950:2016
  LCfitMale <- fit(LC, data = EWStMoMoMale, ages.fit = ages.fit, years.fit = years.fit)
  
  t <- 50
  LCforecastMale <- forecast(LCfitMale, h = t)
  
  mortality_model_age_effect <- LCforecastMale
  
  
  # Since the simulation based on the rates, we add here a row effect
  age_effects <- rep((1+(1:nrow(mortality_model_age_effect$rates))*0.01),
                     ncol(mortality_model_age_effect$rates))
  
  age_effects <- matrix(age_effects, nrow=nrow(mortality_model_age_effect$rates), byrow = FALSE)
  
  mortality_model_age_effect$rates <- age_effects* mortality_model_age_effect$rates
  
  # Super-population
  d_k <- apply(LCforecastMale$rates, 2, function(x) stepfun((min(ages.fit)+1):100, x))
  breaks <- 1:(t-1)
  death_male <- piecewise_xy(breaks,d_k)
  
  # Sub-populations
  d_k_sp <- apply(mortality_model_age_effect$rates, 2, function(x) stepfun((min(ages.fit)+1):100, x))
  death_male_sp <- piecewise_xy(breaks,d_k_sp)
  
  group_effect <- c(1,1.3,0.8)
  
  params <- list("death_male" = death_male, "alpha" = group_effect)
  params_sp <- list("death_male" = death_male_sp, "alpha" = group_effect)
  
  params$mu <- c(0.02,0.001,0.06)
  params_sp$mu <- c(0.02,0.001,0.06)
  
  params$lambda <- 10000 # Entry events
  params_sp$lambda <- 10000
  
  
  death_event <- mk_event_individual(
    type = "death",
    intensity_code = "result = alpha[I.risk_cls-1] * death_male(t,age(I, t));"
  )
  
  exit_event <- mk_event_individual(
    type = "exit",
    intensity = "result = mu[I.risk_cls-1]; "
  )
  
  params$p <- params_sp$p <- 0.5
  
  
  
  entry_event <- mk_event_poisson(
    type = "entry",
    intensity = "lambda",
    kernel_code = "if (CUnif()<p)
newI.risk_cls =2;
else
newI.risk_cls= 3;
double a = CUnif(65,70);
newI.set_age(a,t);
newI.out = false;"
  )
  
  
  model <- mk_model(
    characteristics = get_characteristics(superpop_df), # Characteristics names and types
    events = list(death_event,entry_event, exit_event), # Events list
    parameters = params # Model parameters
  )
  
  model_sp <- mk_model(
    characteristics = get_characteristics(subpop_df), # Characteristics names and types
    events = list(death_event,entry_event, exit_event), # Events list
    parameters = params_sp # Model parameters
  )
  
  death_max <- death_max_sp <- max(sapply(d_k, function(x) { max(x) }))
  
  {
    set.seed(seed)
    sim_out_superp <- popsim(model = model,
                      initial_population = superpop_df,
                      events_bounds = c('death' = death_max, "entry"=params$lambda, "exit"=max(params$mu)),
                      parameters = params,
                      time = t,
                      age_max = 110,
                      multithreading = TRUE)
    
    sim_out_subp <- popsim(model = model_sp,
                             initial_population = subpop_df,
                             events_bounds = c('death' = death_max_sp, "entry"=params_sp$lambda, "exit"=max(params_sp$mu)),
                             parameters = params_sp,
                             time = t,
                             age_max = 110,
                             multithreading = TRUE)
    
    }
  
  
  sim_out <- list(sim_out_superp=sim_out_superp,
                  sim_out_subp=sim_out_subp
                  )
  
  return(sim_out)
  
  
}



# Group and age effect -----


data_generator_1 <- function(age_sample_size = 100000,
                             super_pop_share = 0.8,
                             seed = 1996){
  
  N <- age_sample_size
  
  
  Q1 <- super_pop_share
  N1 <- round(N*Q1)
  ages_1 <- -80:0
  tmp_df <- data.frame("birth"=rep(ages_1,each=N),
                       "death"=rep(rep(NA,N),each=length(ages_1)),
                       "risk_cls"= 1)
  
  ages <- -65:-50
  age_classes <- length(ages)
  
  pop_df <- data.frame("birth"=rep(ages,each=N-N1),
                       "death"=rep(rep(NA,N-N1),each=age_classes),
                       "risk_cls"= sample(c(2,3), 
                                          size=(N-N1)*age_classes,
                                          replace=T,
                                          prob=c(.8,.2)))
  
  pop_df <- rbind(tmp_df,
                  pop_df)
  
  
  
  
  pop_init <- population(pop_df, entry=TRUE, out=TRUE)
  
  EWStMoMoMale <- StMoMoData(EWdata_hmd, series = "male")
  
  #Fitting
  LC <- lc()
  ages.fit <- 50:100
  years.fit <- 1950:2016
  LCfitMale <- fit(LC, data = EWStMoMoMale, ages.fit = ages.fit, years.fit = years.fit)
  
  t <- 50
  LCforecastMale <- forecast(LCfitMale, h = t)
  plot(LCforecastMale)
  
  
  d_k <- apply(LCforecastMale$rates, 2, function(x) stepfun((min(ages.fit)+1):100, x))
  breaks <- 1:(t-1)
  death_male <- piecewise_xy(breaks,d_k)
  
  params <- list("death_male" = death_male, "alpha" = c(1.3,0.8,1))
  
  params$mu <- c(0.001,0.06,0.02)
  
  params$lambda <- 10000 # Entry events
  
  # params$lambda <- 300000 # Entry events
  
  
  death_event <- mk_event_individual(
    type = "death",
    intensity_code = "result = alpha[I.risk_cls-1] * death_male(t,age(I, t));"
  )
  
  exit_event <- mk_event_individual(
    type = "exit",
    intensity = "result = mu[I.risk_cls-1]; "
  )
  
  params$p <- 0.5
  
  
  entry_event <- mk_event_poisson(
    type = "entry",
    intensity = "lambda",
    kernel_code = "if (CUnif()<p)
newI.risk_cls =2;
else
newI.risk_cls= 3;
double a = CUnif(65,70);
newI.set_age(a,t);
newI.out = false;"
  )
  
  model <- mk_model(
    characteristics = get_characteristics(pop_init), # Characteristics names and types
    events = list(death_event,entry_event, exit_event), # Events list
    parameters = params # Model parameters
  )
  
  death_max <- max(sapply(d_k, function(x) { max(x) }))
  
  {
    set.seed(seed)
    sim_out <- popsim(model = model,
                      initial_population = pop_init,
                      events_bounds = c('death' = death_max, "entry"=params$lambda, "exit"=max(params$mu)),
                      parameters = params,
                      time = t,
                      age_max = 110,
                      multithreading = TRUE)
    
  }
  
  
  return(sim_out)
  
  
}

