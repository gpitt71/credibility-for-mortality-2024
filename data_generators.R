# Libraries  ----
rm(list=ls())
Sys.setenv(LANGUAGE = "en")
library(IBMPopSim)
library(dplyr)
library(StMoMo)
library(ggplot2)
library(tidyr)
library(data.table)

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


data_generator_1 <- function(super_pop_cohort_size = 100000,
                             sub_pops_cohort_size=10000,
                             seed = 1996,
                             risk_classes_distribution=c(.6,.4)){
  
  
  
  N1 <- super_pop_cohort_size
   
  ages_1 <- -80:0
  tmp_df <- data.frame("birth"=rep(ages_1,each=N1),
                       "death"=rep(rep(NA,N1),each=length(ages_1)),
                       "risk_cls"= 1)
  
  ages <- (-60:-10)
  age_classes <- length(ages)
  births=rep(ages,each=sub_pops_cohort_size)
  
  
  
  
  pop_df <- data.frame("birth"=births,
                       "death"=NA,#rep(rep(NA,N-N1),each=age_classes),
                       "risk_cls"= sample(c(2,3), 
                                          size=sub_pops_cohort_size*age_classes,
                                          replace=T,
                                          prob=risk_classes_distribution))
  
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
  d_k <- apply(LCforecastMale$rates, 2, function(x) stepfun((min(ages.fit)+1):100, x))
  breaks <- 1:(t-1)
  death_male <- piecewise_xy(breaks,d_k)
  
  params <- list("death_male" = death_male, "alpha" = c(1,0.8,1.3))
  
  params$mu <- c(0.001,0.006,0.002)
  
  params$lambda <- 1000 # Entry events
  
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
double a = CUnif(20,60);
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



data_generator_hmd <- function(seed_input=1,
                               exposure_sup=100000,
                               exposure_sub=c(5000,
                                              500)){
  
  
  
  
  dt <- read.table(
    "C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\data\\lt_ITA.txt",
    header = T
  ) %>%
    rename(Cohort = Year) %>%
    select(Cohort, Age, qx,lx) %>%
    mutate(
      Age = replace(Age, Age == "110+", "110"),
      Age = as.integer(Age),
      qx = as.numeric(qx))%>%
    mutate_all(~replace(., is.na(.) | is.infinite(.), 1)) %>%
  mutate(intercept = log(qx / (1 - qx))) %>%
    mutate_all(~replace(., is.infinite(.), 0))
  
  
  
  {
    set.seed(seed = seed_input)
    age_eff_1 <- 1 - runif(length(0:110), .2, .3)
    age_eff_2 <- 1 + runif(length(0:110), .2, .3)
    }
  
  dt_theta <- data.frame(Age = 0:110,
                         Theta1x = age_eff_1,
                         Theta2x = age_eff_2)
  
  
  dt <- dt %>%
    mutate(Cohort = Cohort -min(Cohort) , Period = Cohort + Age) %>%
    filter(Cohort >= 0)
  
  dt <- left_join(dt, dt_theta, by = "Age") %>%
    mutate(
      predictor1 = Theta1x * exp(intercept),
      predictor2 = Theta2x * exp(intercept),
      p1 = (predictor1 / (1 + predictor1)),
      p2 = (predictor2 / (1 + predictor2)),
      p = qx
    )
  
  
  
  
  dt <- dt%>%
    select(Age,Period,Cohort,p,p1,p2) 
  
  
  setDT(dt)  
  
  
  dt <- dt %>%
    arrange(Cohort, Age, Period)
  
  dt[,c('Lx','Dx'):=binomial_simulator(starting_exposure=exposure_sup,
                                       probabilities = p),
     by=.(Cohort)]  
  
  
  dt[,c('Lx1','Dx1'):=binomial_simulator(starting_exposure=exposure_sub[1],
                                       probabilities = p1),
     by=.(Cohort)]  
  
  
  dt[,c('Lx2','Dx2'):=binomial_simulator(starting_exposure=exposure_sub[2],
                                         probabilities = p2),
     by=.(Cohort)]  
  
  dt <- dt %>%
    arrange(Cohort, Age, Period)
  
  
  dt[,c('Ext','Dx'):=average_exposure_and_deaths(raw_exposure=Lx),
     by=.(Cohort)]  
  
  dt[,c('Ext1','Dx1'):=average_exposure_and_deaths(raw_exposure=Lx1),
     by=.(Cohort)]  
  
  dt[,c('Ext2','Dx2'):=average_exposure_and_deaths(raw_exposure=Lx2),
     by=.(Cohort)]  
  
  
  return(dt)
  
}




data_generator_hmd_2 <- function(seed_input=1,
                               exposure_sup=100000,
                               exposure_sub=c(5000,
                                              500)){
  
  
  dt <- read.table(
    "C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\data\\lt_ITA.txt",
    header = T
  ) %>%
    rename(Cohort = Year) %>%
    select(Cohort, Age, qx,lx) %>%
    mutate(
      Age = replace(Age, Age == "110+", "110"),
      Age = as.integer(Age),
      qx = as.numeric(qx))%>%
    mutate_all(~replace(., is.na(.) | is.infinite(.), 1)) %>%
    mutate(intercept = log(qx / (1 - qx))) %>%
    mutate_all(~replace(., is.infinite(.), 0))
  
  
  
  {
    set.seed(seed = seed_input)
    Theta1x=1 - runif(1, .2, .3)
    Theta2x=1 + runif(1, .2, .3)
    # age_eff_1 <- 1 - runif(length(0:110), .2, .3)
    # age_eff_2 <- 1 + runif(length(0:110), .2, .3)
    }
  
  # dt_theta <- data.frame(Age = 0:110,
  #                        Theta1x = age_eff_1,
  #                        Theta2x = age_eff_2)
  
  
  dt <- dt %>%
    mutate(Cohort = Cohort -min(Cohort) , Period = Cohort + Age) %>%
    filter(Cohort >= 0)
  
  dt <- dt%>%#left_join(dt, dt_theta, by = "Age") %>%
    mutate(
      predictor1 = Theta1x * exp(intercept),
      predictor2 = Theta2x * exp(intercept),
      p1 = (predictor1 / (1 + predictor1)),
      p2 = (predictor2 / (1 + predictor2)),
      p = qx
    )
  
  
  
  
  dt <- dt%>%
    select(Age,Period,Cohort,p,p1,p2) 
  
  
  setDT(dt)  
  
  
  dt <- dt %>%
    arrange(Cohort, Age, Period)
  
  dt[,c('Lx','Dx'):=binomial_simulator(starting_exposure=exposure_sup,
                                       probabilities = p),
     by=.(Cohort)]  
  
  
  dt[,c('Lx1','Dx1'):=binomial_simulator(starting_exposure=exposure_sub[1],
                                         probabilities = p1),
     by=.(Cohort)]  
  
  
  dt[,c('Lx2','Dx2'):=binomial_simulator(starting_exposure=exposure_sub[2],
                                         probabilities = p2),
     by=.(Cohort)]  
  
  dt <- dt %>%
    arrange(Cohort, Age, Period)
  
  
  dt[,c('Ext','Dx'):=average_exposure_and_deaths(raw_exposure=Lx),
     by=.(Cohort)]  
  
  dt[,c('Ext1','Dx1'):=average_exposure_and_deaths(raw_exposure=Lx1),
     by=.(Cohort)]  
  
  dt[,c('Ext2','Dx2'):=average_exposure_and_deaths(raw_exposure=Lx2),
     by=.(Cohort)]  
  
  
  return(dt)
  
}




data_generator_hmd_completing_the_cohorts <- function(seed_input=1,
                               exposure_sup=100000,
                               exposure_sub=c(5000,
                                              500)){
  
  
  
  
  dt <- read.table(
    "C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\data\\lt_ITA.txt",
    header = T
  ) %>%
    rename(Cohort = Year) %>%
    select(Cohort, Age, qx,lx) %>%
    mutate(
      Age = replace(Age, Age == "110+", "110"),
      Age = as.integer(Age),
      qx = as.numeric(qx))%>%
    mutate_all(~replace(., is.na(.) | is.infinite(.), 1)) %>%
    mutate(intercept = log(qx / (1 - qx))) %>%
    mutate_all(~replace(., is.infinite(.), 0))
  
  
  
  {
    set.seed(seed = seed_input)
    age_eff_1 <- 1 - runif(length(0:110), .2, .3)
    age_eff_2 <- 1 + runif(length(0:110), .2, .3)
    }
  
  dt_theta <- data.frame(Age = 0:110,
                         Theta1x = age_eff_1,
                         Theta2x = age_eff_2)
  
  
  dt <- dt %>%
    mutate(Cohort = Cohort -min(Cohort) , Period = Cohort + Age) %>%
    filter(Cohort >= 0)
  
  dt <- left_join(dt, dt_theta, by = "Age") %>%
    mutate(
      predictor1 = Theta1x * exp(intercept),
      predictor2 = Theta2x * exp(intercept),
      p1 = (predictor1 / (1 + predictor1)),
      p2 = (predictor2 / (1 + predictor2)),
      p = qx
    )
  
  
  
  
  dt <- dt%>%
    select(Age,Period,Cohort,p,p1,p2) 
  
  
  setDT(dt)  
  
  
  dt <- dt %>%
    arrange(Cohort, Age, Period)
  
  dt[,c('Lx','Dx'):=binomial_simulator(starting_exposure=exposure_sup,
                                       probabilities = p),
     by=.(Cohort)]  
  

  
  dt[,c('Ext','Dx'):=average_exposure_and_deaths(raw_exposure=Lx),
     by=.(Cohort)]  
  
  tmp_occ <-  dt %>%
    tidyr::pivot_wider(id_cols=Age,names_from = Period, values_from = Dx) %>%
    as.matrix()
  rownames(tmp_occ) <- tmp_occ[,"Age"]
  tmp_occ <- tmp_occ[,-1]
  # tmp_occ <- tmp_occ[as.character(ages_fit),years_fit+bias]
  
  
  tmp_exp <-  dt %>%
    tidyr::pivot_wider(id_cols=Age,names_from = Period, values_from = Ext) %>%
    as.matrix()
  rownames(tmp_exp) <- tmp_exp[,"Age"]
  tmp_exp <- tmp_exp[,-1]
  # tmp_exp <- tmp_exp[as.character(ages_fit),years_fit+bias]
  # colnames(tmp_occ) <-colnames(tmp_exp) <- years_fit-min(years_fit)
  
  
  datahat <- structure(
    list(
      Dxt = tmp_occ,
      Ext = tmp_exp,
      ages = as.numeric(rownames(tmp_exp)),
      years = as.numeric(colnames(tmp_exp)),
      type = 'central',
      series = 'female',
      label = 'total'
    ),
    class = "StMoMoData"
  )
  
  ap = StMoMo::StMoMo(link="logit",
                      staticAgeFun = TRUE,
                      periodAgeFun = c("1"),
                      cohortAgeFun = NULL)
  
  model_fit <- fit(ap, data = central2initial(datahat))
  
  
  
  cal_component <- matrix(rep(model_fit$kt, nrow(tmp_exp)),
                          nrow = nrow(tmp_exp),
                          byrow = T)
  
  age_component <- matrix(rep(model_fit$ax, ncol(tmp_exp)),
                          nrow = nrow(tmp_exp),
                          byrow = F)
  
  
  mx <- exp(age_component+cal_component)/(1+exp(age_component+cal_component))
  
  mx <- data.frame(mx)
  
  colnames(mx) <- datahat$years
  mx[,'Age'] <- datahat$ages
  
  mx <- pivot_longer(mx,
                     cols=as.character(0:168),
                     names_to = "Period",
                     values_to = "p"
                     )
  
  mx <- mx %>%
    mutate(Period=as.integer(Period),
           Cohort= Period-Age)
  
  missing_cohorts <- setdiff(mx$Cohort,dt$Cohort)
  
  mx <- mx %>% filter(Cohort %in% missing_cohorts)
  
  dt <- rbind(dt %>%
                select(Age, Period, Cohort, p), mx)
  
  setDT(dt)  
  
  dt <- left_join(dt, dt_theta, by = "Age") %>%
    mutate(
      intercept=log(p/(1-p)),
      predictor1 = Theta1x * exp(intercept),
      predictor2 = Theta2x * exp(intercept),
      p1 = (predictor1 / (1 + predictor1)),
      p2 = (predictor2 / (1 + predictor2))
    )
  
  
  dt <- dt %>%
    arrange(Cohort, Age, Period)
  
  dt[,c('Lx','Dx'):=binomial_simulator(starting_exposure=exposure_sup,
                                       probabilities = p),
     by=.(Cohort)]  
  
  
  dt[,c('Lx1','Dx1'):=binomial_simulator(starting_exposure=exposure_sub[1],
                                         probabilities = p1),
     by=.(Cohort)]
  
  
  dt[,c('Lx2','Dx2'):=binomial_simulator(starting_exposure=exposure_sub[2],
                                         probabilities = p2),
     by=.(Cohort)]
  
  dt <- dt %>%
    arrange(Cohort, Age, Period)
  

  
  dt[,c('Ext','Dx'):=average_exposure_and_deaths(raw_exposure=Lx),
     by=.(Cohort)]  
  
  dt[,c('Ext1','Dx1'):=average_exposure_and_deaths(raw_exposure=Lx1),
     by=.(Cohort)]  
  
  dt[,c('Ext2','Dx2'):=average_exposure_and_deaths(raw_exposure=Lx2),
     by=.(Cohort)]  
  
  
  return(dt)
  
}





data_generator_hmd_periods <- function(seed_input=1,
                                                      exposure_sup=100000,
                                                      exposure_sub=c(5000,
                                                                     500)){
  
  
  
  
  dt <- read.table(
    "C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\data\\lt_ITA_period.txt",
    header = T
  ) %>%
    rename(Period = Year) %>%
    select(Period, Age, qx,lx) %>%
    mutate(
      Age = replace(Age, Age == "110+", "110"),
      Age = as.integer(Age),
      qx = as.numeric(qx))%>%
    mutate_all(~replace(., is.na(.) | is.infinite(.), 1)) %>%
    mutate(intercept = log(qx / (1 - qx))) %>%
    mutate_all(~replace(., is.infinite(.), 0))
  
  
  
  {
    set.seed(seed = seed_input)
    age_eff_1 <- 1 - runif(length(0:110), .2, .3)
    age_eff_2 <- 1 + runif(length(0:110), .2, .3)
    }
  
  dt_theta <- data.frame(Age = 0:110,
                         Theta1x = age_eff_1,
                         Theta2x = age_eff_2)
  
  
  dt <- dt %>%
    mutate(Period = Period -min(Period) , Cohort = Period-Age)%>%
    filter(Cohort >= 0)

  dt <- left_join(dt, dt_theta, by = "Age") %>%
    mutate(
      predictor1 = Theta1x * exp(intercept),
      predictor2 = Theta2x * exp(intercept),
      p1 = (predictor1 / (1 + predictor1)),
      p2 = (predictor2 / (1 + predictor2)),
      p = qx
    )

  
  dt <- dt%>%
    select(Age,Period,Cohort,p,p1,p2) 
  
  
  setDT(dt)  
  
  
  dt <- dt %>%
    arrange(Cohort, Age, Period)
  
  dt[,c('Lx','Dx'):=binomial_simulator(starting_exposure=exposure_sup,
                                       probabilities = p),
     by=.(Cohort)]  
  
  dt[,Dx:=as.numeric(Dx)]
  
  
  dt[,c('Lx1','Dx1'):=binomial_simulator(starting_exposure=exposure_sub[1],
                                         probabilities = p1),
     by=.(Cohort)]  
  
  dt[,Dx1:=as.numeric(Dx1)]
  
  dt[,c('Lx2','Dx2'):=binomial_simulator(starting_exposure=exposure_sub[2],
                                         probabilities = p2),
     by=.(Cohort)]  
  
  dt[,Dx2:=as.numeric(Dx2)]
  
    dt <- dt %>%
    arrange(Cohort, Age, Period)
  
  
  dt[,c('Ext','Dx'):=average_exposure_and_deaths(raw_exposure=Lx),
     by=.(Cohort)]  
  
  dt[,c('Ext1','Dx1'):=average_exposure_and_deaths(raw_exposure=Lx1),
     by=.(Cohort)]  
  
  dt[,c('Ext2','Dx2'):=average_exposure_and_deaths(raw_exposure=Lx2),
     by=.(Cohort)]  
  
  
  return(dt)
  
}


data_generator_hmd_lt <- function(seed_input=1,
                                  exposure_sup=100000,
                                  exposure_sub=c(5000,
                                                 500)){
  
  
  ita_data <- demography::hmd.mx(country = "ITA",
                  "gabriele.pittarello@uniroma1.it",
                  "pm869fqW8nxozS8!")
  
  
  E1 <- ita_data$pop[['female']]

  D1 <- E1  * ita_data$rate[['female']]
  
  periods_cols <- as.numeric(colnames(E1))
  
  
  periods_cols <- periods_cols - min(periods_cols)
  
  colnames(E1) <- colnames(D1) <- as.character(periods_cols)
  
  datahat <- structure(
    list(
      Dxt = D1,
      Ext = E1,
      ages = as.numeric(rownames(E1)),
      years = as.numeric(colnames(E1)),
      type = 'central',
      series = 'female',
      label = 'total'
    ),
    class = "StMoMoData"
  )
  
  initial_data <- central2initial(datahat)
  
  tmp_occ <- initial_data$Dxt
  tmp_exp <- initial_data$Ext
  
  
  tmp_occ <- as.data.frame(tmp_occ)
  tmp_exp <- as.data.frame(tmp_exp)
  tmp_exp[['Age']] <-tmp_occ[['Age']] <- as.numeric(rownames(tmp_occ))
  
  dt1 <- pivot_longer(tmp_occ,
                      cols = as.character(periods_cols),
                      names_to="Period",
                      values_to = "Dxt")

  dt2 <- pivot_longer(tmp_exp,
                      cols = as.character(periods_cols),
                      names_to="Period",
                      values_to = "Ext")
  

  dt <- left_join(dt1,dt2,by=c("Age","Period"))
  
  dt <- dt %>% mutate(qx = Dxt/Ext,
                      Period=as.numeric(Period),
                      Cohort = Period-Age) 
  
  
  # dt <- dt %>% 
    # select(Age,Period,Cohort,qx) %>%
    # fill_missing_rows()
  
  dt <- dt%>%
    mutate_all(~replace(., is.na(.) | is.infinite(.), 1)) %>%
    mutate(intercept = log(qx / (1 - qx))) %>%
    mutate_all(~replace(., is.infinite(.), 0))
  
  

  
  {
    set.seed(seed = seed_input)
    age_eff_1 <- 1 - runif(length(0:110), .2, .3)
    age_eff_2 <- 1 + runif(length(0:110), .2, .3)
    }
  
  dt_theta <- data.frame(Age = 0:110,
                         Theta1x = age_eff_1,
                         Theta2x = age_eff_2)
  
  
  # dt <- dt %>%
  #   mutate(Period = Period -min(Period) , Cohort = Period-Age)%>%
    # filter(Cohort >= 0)
  # 
  
  dt <- left_join(dt, dt_theta, by = "Age") %>%
    mutate(
      predictor1 = Theta1x * exp(intercept),
      predictor2 = Theta2x * exp(intercept),
      p1 = (predictor1 / (1 + predictor1)),
      p2 = (predictor2 / (1 + predictor2)),
      p = qx
    )
  
  
  dt <- dt%>%
    select(Age,Period,Cohort,p,p1,p2) 
  
  
  setDT(dt)  
  
  
  dt <- dt %>%
    arrange(Cohort, Age, Period)
  
  dt[,c('Lx','Dx_tmp'):=binomial_simulator(starting_exposure=exposure_sup,
                                       probabilities = p),
     by=.(Cohort)]  
  
  dt[,Dx_tmp:=as.numeric(Dx_tmp)]
  
  
  dt[,c('Lx1','Dx1_tmp'):=binomial_simulator(starting_exposure=exposure_sub[1],
                                         probabilities = p1),
     by=.(Cohort)]  
  
  dt[,Dx1_tmp:=as.numeric(Dx1_tmp)]
  
  dt[,c('Lx2','Dx2_tmp'):=binomial_simulator(starting_exposure=exposure_sub[2],
                                         probabilities = p2),
     by=.(Cohort)]  
  
  dt[,Dx2_tmp:=as.numeric(Dx2_tmp)]
  
  dt <- dt %>%
    arrange(Cohort, Age, Period)
  
  
  dt[,c('Ext','Dx'):=average_exposure_and_deaths(raw_exposure=Lx),
     by=.(Cohort)]  
  
  dt[,c('Ext1','Dx1'):=average_exposure_and_deaths(raw_exposure=Lx1),
     by=.(Cohort)]  
  
  dt[,c('Ext2','Dx2'):=average_exposure_and_deaths(raw_exposure=Lx2),
     by=.(Cohort)]  
  
  
  return(dt)
  
  
  
  }










