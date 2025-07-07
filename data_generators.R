# Libraries  ----
rm(list=ls())
Sys.setenv(LANGUAGE = "en")
library(dplyr)
library(StMoMo)
library(ggplot2)
library(tidyr)
library(data.table)
library(rpart)

# Group and age effect -----

data_generator_hmd_lt <- function(seed_input=1,
                                  hmd_username=NULL,
                                  hmd_password=NULL,
                                  exposure_sup=100000,
                                  exposure_sub=c(5000,
                                                 500)){
  
  
  ita_data <- demography::hmd.mx(country = "ITA",
                                 hmd_username,
                                 hmd_password)
  
  
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


# Print group and age effect ----

print_age_effects <- function(seed_input=1,
                                  hmd_username=NULL,
                                  hmd_password=NULL,
                                  exposure_sup=100000,
                                  exposure_sub=c(5000,
                                                 500)){
  
  
  ita_data <- demography::hmd.mx(country = "ITA",
                                 hmd_username,
                                 hmd_password)
  
  
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
  
  
  out<-list(age_eff_1=age_eff_1,
       age_eff_2=age_eff_2)
 
  return(out)
   
}







