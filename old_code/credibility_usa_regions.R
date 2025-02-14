rm(list = ls())
library(StMoMo)
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(purrr)
library(demography)

# Read data

files <- stringr::str_glue(
  "C:\\Users\\gpitt\\Documents\\Postdoc\\Torino\\Mortality\\data\\usa_divisions\\Div{1:9}_usa.csv"
)

dt <-  files %>%
  map_dfr(readr::read_csv)

dt = dt %>%
  mutate(Ext = dx / mx)

list_of_dfs <- dt %>% split(dt$PopName)

create_df_usa_division <- function(x) {
  exposure <- pivot_wider(x,
                          id_cols = "Age",
                          values_from = "Ext",
                          names_from = "Year")
  
  exposure <- exposure %>% tibble::column_to_rownames("Age")
  
  
  
  rates <- pivot_wider(x,
                       id_cols = "Age",
                       values_from = "mx",
                       names_from = "Year")
  
  rates <- rates %>% tibble::column_to_rownames("Age")
  
  deaths <- pivot_wider(x,
                       id_cols = "Age",
                       values_from = "dx",
                       names_from = "Year")
  
  deaths <- deaths %>% tibble::column_to_rownames("Age")
  
  l <- list(exposure = as.matrix.data.frame(exposure), 
            rates = as.matrix.data.frame(rates),
            deaths=as.matrix.data.frame(deaths))
  
  return(l)
  
  
}


l <- lapply(list_of_dfs, 
            create_df_usa_division)

ages.fit <- 55:89
years_fit_basic <- 1951:2000

series_code='female'


mse_0 <- mse_1 <- mse_2 <- NULL


prediction_horizon <- 0:18

for (ix in prediction_horizon) {
  years.fit <- min(years_fit_basic):max(years_fit_basic + ix)
  
  
  for (i in 1:9) {
    assign(paste0("E", i), l[[i]]$exposure[as.character(ages.fit), as.character(years.fit)])
    assign(paste0("muxt", i), l[[i]]$rates[as.character(ages.fit), as.character(years.fit)])
    assign(paste0("D", i), l[[i]]$deaths[as.character(ages.fit), as.character(years.fit)])
  }
  
  Ext <- 0
  Dxt <- 0
  
  for (i in 1:9) {
    Ext <- Ext + get(paste0("E", i))
    Dxt <- Dxt + get(paste0("D", i))
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
  # LC <- lc(link = "log")
  
  LC <- apc(link="log")

    LCModelfit <- fit(LC,
                      data = datahat,
                      years.fit = years.fit,
                      ages.fit = ages.fit)

    year.predict <- max(years.fit)+1

    LCModelfor <- forecast(LCModelfit,
                           h = 1)
    
    
    for (i in 1:9) {
      assign(paste0("muxt_actual_", i), l[[i]]$rates[as.character(ages.fit), as.character(year.predict)])
      
      }
  

      muxt_hat_predicted <- LCModelfor$rates
      
      
      mse_0_ix <- 0
      
      for (i in 1:9) {
        mse_0_ix <- mse_0_ix + (get(paste0("muxt_actual_", i))-muxt_hat_predicted)^2
      }
      
      mse_0 <- c(mse_0,sum(mse_0_ix))
      
      
      muxt_hat <- LCModelfor$fitted
      

      for (i in 1:9) {
        assign(paste0("Fxt_", i), get(paste0("D",i))/get(paste0("E",i)))
        assign(paste0("varthetax_", i), apply(((get(paste0("Fxt_",i))-muxt_hat)^2)/(muxt_hat^2),1,mean))
        assign(paste0("C", i), apply(get(paste0("D",i)),1,sum)/(apply(get(paste0("E",i))*muxt_hat,1,sum)))
        assign(paste0("Z_", i), 1/(1+get(paste0("varthetax_", i))*apply(get(paste0("E",i))*muxt_hat,1,sum)))
        assign(paste0("muhat", i), get(paste0("Z_", i))*muxt_hat_predicted+(1-get(paste0("Z_", i)))*get(paste0("C", i))*muxt_hat_predicted)
      }
      
      mse_1_ix <- 0
      
      for (i in 1:9) {
        mse_1_ix <- mse_1_ix + (get(paste0("muxt_actual_", i))-get(paste0("muhat", i)))^2
      }
      
      mse_1 <- c(mse_1,sum(mse_1_ix))
      
      mse_2_ix <- 0
      
      for(i in 1:9){
        
          popspecdata <- structure(
            list(
              Dxt = l[[i]]$deaths[as.character(ages.fit), as.character(years.fit)],
              Ext = l[[i]]$exposure[as.character(ages.fit), as.character(years.fit)],
              ages = ages.fit,
              years = years.fit,
              type = 'central',
              series = series_code,
              label =  paste0("Div",i)
            ),
            class = "StMoMoData"
          )
          
          assign(paste0("LCModelfit",i),
                 fit(LC,
                     data = popspecdata,
                     years.fit = years.fit,
                     ages.fit = ages.fit))
          
          assign(paste0("LCModelfor",i),forecast(get(paste0("LCModelfit",i)), h = 1))
          
          mse_2_ix <- mse_2_ix + (get(paste0("muxt_actual_", i))-get(paste0("LCModelfor", i))$rates)^2
          
      }
      
      mse_2 <- c(mse_2,sum(mse_2_ix))
    
}


print(c(sum(mse_0)/(length(ages.fit)*9*length(prediction_horizon))*10e+03,
      sum(mse_1)/(length(ages.fit)*9*length(prediction_horizon))*10e+03,
      sum(mse_2)/(length(ages.fit)*9*length(prediction_horizon))*10e+03))
