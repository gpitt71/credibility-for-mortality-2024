rm(list = ls())
library(StMoMo)
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(purrr)
library(demography)
library(data.table)

# Read data


tmp <- fread("C:\\Users\\gpitt\\Documents\\Postdoc\\Torino\\Mortality\\data\\ita_regions\\region1.csv")


tmp <- tmp %>%
  rename("Age"="Età e classi di età",
         "Period"="TIME",
         "Region"="Territorio",
         "Gender"="Sesso",
         "Information"="Funzioni biometriche") %>%
  select(Age,Period,Region,Gender,Information,Value)%>%
  arrange(Period)

# mutate(Region = case_when(
#   Region == "Piemonte" ~ "Piemonte_VDA",   
#   Region == "Valle d'Aosta / Vallée d'Aoste" ~ "Piemonte_VDA",   
#   Region == "Molise" ~ "Marche_Molise",   
#   Region == "Marche" ~ "Marche_Molise",   
#   TRUE ~ Region  
# )) %>%

tmp <- tmp %>%
  filter(Region != "Valle d'Aosta / Vallée d'Aoste" & 
           Region != "Molise" )
           
           
tmp <- tmp %>%
  mutate(Age=as.integer(gsub("([0-9]+).*$", "\\1", Age)))
  
long_format <- pivot_wider(tmp,
                           id_cols = c("Age",
                                       "Period",
                                       "Region",
                                       "Gender"),
                           names_from="Information",
                           values_from = "Value")




long_format <- long_format %>% 
  filter(Gender=="femmine")

colnames(long_format)[5] <- "lx"
colnames(long_format)[6] <- "dx"

long_format <- long_format %>% 
  arrange(Period,Age) %>%
  group_by(Region,Period) %>%
  mutate(lx_lead=lead(lx))

long_format <- long_format %>% 
  as.data.frame() %>%
  mutate(mx = dx/((lx+lx_lead)/2),
         Ext = (lx+lx_lead)/2)

list_of_dfs <- long_format %>% split(long_format$Region)



create_df_ita_regions <- function(x) {
  exposure <- pivot_wider(x,
                          id_cols = "Age",
                          values_from = "Ext",
                          names_from = "Period")
  
  exposure <- exposure %>% tibble::column_to_rownames("Age")
  
  
  
  rates <- pivot_wider(x,
                       id_cols = "Age",
                       values_from = "mx",
                       names_from = "Period")
  
  rates <- rates %>% tibble::column_to_rownames("Age")
  
  deaths <- pivot_wider(x,
                        id_cols = "Age",
                        values_from = "dx",
                        names_from = "Period")
  
  deaths <- deaths %>% tibble::column_to_rownames("Age")
  
  l <- list(exposure = as.matrix.data.frame(exposure), 
            rates = as.matrix.data.frame(rates),
            deaths=as.matrix.data.frame(deaths))
  
  return(l)
  
  
}


l <- lapply(list_of_dfs, 
            create_df_ita_regions)

number_of_groups <- length(l)
N_groups <- length(l)

ages.fit <- 50:89
years_fit_basic <- 1974:2003

series_code='male'

mortality_model_lc <- lc(link="log")
mortality_model_apc <- apc(link="log")
mortality_model_rh <- rh(link="log")

out <- NULL 
for(forecasting_horizon in c(2,8,12)){
  
  prediction_horizon <- 0:(15-forecasting_horizon)
  val_years = forecasting_horizon-1

  for(model_option in c("lc","apc","rh")){
    
    mse_0_val <- mse_1_val <- mse_2_val <- NULL
    mae_0_val <- mae_1_val <- mae_2_val <- NULL
    
    mse_0 <- mse_1 <- mse_2 <- NULL
    mae_0 <- mae_1 <- mae_2 <- NULL
    
    assign("mortality_model",get(paste0("mortality_model_",model_option)))
    

for (ix in prediction_horizon) {
  years.fit <- min(years_fit_basic):max(years_fit_basic + ix)
  
  
  for (i in 1:number_of_groups) {
    assign(paste0("E", i), l[[i]]$exposure[as.character(ages.fit), as.character(years.fit)])
    assign(paste0("muxt", i), l[[i]]$rates[as.character(ages.fit), as.character(years.fit)])
    assign(paste0("D", i), l[[i]]$deaths[as.character(ages.fit), as.character(years.fit)])
  }
  
  Ext <- 0
  Dxt <- 0
  
  # for (i in 1:number_of_groups) {
  #   Ext <- Ext + get(paste0("E", i))
  #   Dxt <- Dxt + get(paste0("D", i))
  # }
  
  # for (i in 1:number_of_groups) {
  #   Ext <- Ext + get(paste0("E", i))
  #   Dxt <- Dxt + get(paste0("D", i))
  # }
  
  list_of_extra_exposures <- list()
  
  for (i in 2:number_of_groups) {
    
    tmp_list <- list()
    
    tmp_list[['Dxt']] <- get(paste0("D", i))
    tmp_list[['Ext']] <- get(paste0("E", i))
    list_of_extra_exposures[[i-1]] <- tmp_list
    
  }

  
  datahat <- structure(
    list(
      Dxt = D1,
      Ext = E1,
      ages = ages.fit,
      years = years.fit,
      type = 'central',
      series = series_code,
      label = 'total'
    ),
    class = "StMoMoData"
  )


  mortality_model_fit <- fit(mortality_model,
                             data = datahat,
                             years.fit = years.fit,
                             ages.fit = ages.fit,
                             list_of_extra_exposures = list_of_extra_exposures)

    year.predict <- max(years.fit)+forecasting_horizon
    
    years.validation <- (max(years.fit)+1):(year.predict-1)
    
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
    
    
    for (i in 1:number_of_groups) {
      assign(paste0("muxt_actual_", i), l[[i]]$rates[as.character(ages.fit), as.character(year.predict)])
      assign(paste0("muxt_actual_", i, "_val"), l[[i]]$rates[as.character(ages.fit), as.character(years.validation)])
      }
  

    muxt_hat_predicted_val <- mortality_model_forecast$rates[as.character(ages.fit),as.character(years.validation)]
    
    
    if(forecasting_horizon>1){
      muxt_hat_predicted <- mortality_model_forecast$rates[as.character(ages.fit),as.character(year.predict)]
    }else{
      muxt_hat_predicted <- mortality_model_forecast$rates
    }
      
      
    mse_0_ix <- 0
    mse_0_ix_val <- 0
    
    for (i in 1:N_groups) {
      mse_0_ix <- mse_0_ix + (abs((get(paste0("muxt_actual_", i))-muxt_hat_predicted)))#/get(paste0("muxt_actual_", i))
      mse_0_ix_val <- mse_0_ix_val + (abs(get(paste0("muxt_actual_", i,"_val"))-muxt_hat_predicted_val))#/get(paste0("muxt_actual_", i,"_val"))
      
    }
    
    mse_0 <- c(mse_0,sum(mse_0_ix))
    mse_0_val <- c(mse_0_val,sum(mse_0_ix_val))
      
      
      muxt_hat <- mortality_model_forecast$fitted
      

      for (i in 1:number_of_groups) {
        assign(paste0("Fxt_", i), get(paste0("D",i))/get(paste0("E",i)))
        assign(paste0("varthetax_", i), apply(((get(paste0("Fxt_",i))-muxt_hat)^2)/(muxt_hat^2),1,mean))
        assign(paste0("C", i), apply(get(paste0("D",i)),1,sum)/(apply(get(paste0("E",i))*muxt_hat,1,sum)))
        assign(paste0("C", i,"_val"), matrix(rep(get(paste0("C",i)),val_years),byrow = F,ncol = val_years))
        
        assign(paste0("Z_", i), 1/(1+get(paste0("varthetax_", i))*apply(get(paste0("E",i))*muxt_hat,1,sum)))
        assign(paste0("Z_", i,"_val"), matrix(rep(get(paste0("Z_",i)),val_years),byrow = F,ncol = val_years))
        
        
        assign(paste0("muhat", i,"_val"),get(paste0("Z_", i,"_val"))*muxt_hat_predicted_val+(1-get(paste0("Z_", i,"_val")))*get(paste0("C", i,"_val"))*muxt_hat_predicted)
        assign(paste0("muhat", i), get(paste0("Z_", i))*muxt_hat_predicted+(1-get(paste0("Z_", i)))*get(paste0("C", i))*muxt_hat_predicted)
      }
      
      mse_1_ix <- 0
      mse_1_ix_val <- 0
      
      for (i in 1:N_groups) {
        
        mse_1_ix_val <- mse_1_ix_val + (abs(get(paste0("muxt_actual_", i,"_val"))-get(paste0("muhat", i,"_val"))))#/get(paste0("muxt_actual_", i,"_val"))
        mse_1_ix <- mse_1_ix + (abs(get(paste0("muxt_actual_", i))-get(paste0("muhat", i))))#/get(paste0("muxt_actual_", i))
        
      }
      
      mse_1 <- c(mse_1,sum(mse_1_ix))
      mse_1_val <- c(mse_1_val,sum(mse_1_ix_val))
      
      mse_2_ix <- 0
      mse_2_ix_val <- 0
      
      for(i in 1:number_of_groups){
        
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
          
          assign(paste0("mortality_model_fit",i),
                 fit(mortality_model,
                     data = popspecdata,
                     years.fit = years.fit,
                     ages.fit = ages.fit))
          
          assign(paste0("cv.arima.kt_", i), auto.arima(as.numeric(get(paste0("mortality_model_fit",i))$kt),
                                                       ic="bic"))
          
          
          if(model_option != "lc"){
            assign(paste0("cv.arima.gc_", i), auto.arima(as.numeric(get(paste0("mortality_model_fit",i))$gc),
                                                         ic="bic"))
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
            assign(paste0("mm_forecast_", i, "_val") , get(paste0("mortality_model_forecast",i))$rates[as.character(ages.fit),as.character(years.validation)])
            
          }else{
            
            assign(paste0("mm_forecast_", i) , get(paste0("mortality_model_forecast",i))$rates)
          }

          mse_2_ix <- mse_2_ix + (abs(get(paste0("muxt_actual_", i))-get(paste0("mm_forecast_", i))))#/get(paste0("muxt_actual_", i))
          mse_2_ix_val <- mse_2_ix_val + (abs(get(paste0("muxt_actual_", i,"_val"))-get(paste0("mm_forecast_", i,"_val"))))#/get(paste0("muxt_actual_", i,"_val"))
          
      }
      
      mse_2 <- c(mse_2,sum(mse_2_ix))
      mse_2_val <- c(mse_2_val,sum(mse_2_ix_val))
    
}
  
    out_model <- data.frame(
      forecasting_horizon = forecasting_horizon,
      model = c(model_option),
      mse_0_val = c(sum(mse_0_val)/(length(ages.fit)*val_years*N_groups*length(prediction_horizon))),
      mse_0 = c(sum(mse_0)/(length(ages.fit)*N_groups*length(prediction_horizon))),
      mse_1_val = c(sum(mse_1_val)/(length(ages.fit)*val_years*N_groups*length(prediction_horizon))),
      mse_1 = sum(mse_1)/(length(ages.fit)*N_groups*length(prediction_horizon)),
      mse_2_val = c(sum(mse_2_val)/(length(ages.fit)*val_years*N_groups*length(prediction_horizon))),
      mse_2 = sum(mse_2)/(length(ages.fit)*N_groups*length(prediction_horizon)))
    
    out <- rbind(out,
                 out_model)
    
    
    
    }}



fwrite(out,
       "C:\\Users\\gpitt\\Documents\\Postdoc\\Torino\\Mortality\\results\\ita_validation_testing_l2.csv")


long_data_val <- pivot_longer(
  out,
  cols = mse_0_val:mse_2,
  names_to = "Metric",
  values_to = "Value"
)%>%
  filter(grepl("_val",Metric)) %>%
  group_by(forecasting_horizon,Metric,model) %>%
  reframe(mse_val=min(Value)) %>%
  group_by(forecasting_horizon,Metric) %>%
  filter(mse_val==min(mse_val)) %>%
  mutate(family = as.numeric(gsub("[^0-9]", "", Metric)))

long_data <- pivot_longer(
  out,
  cols = mse_0_val:mse_2,
  names_to = "Metric",
  values_to = "Value"
) %>%
  filter(!grepl("_val",Metric)) %>%
  mutate(family = as.numeric(gsub("[^0-9]", "", Metric)))



out_selected <-   left_join(long_data_val,
                            long_data,
                            by = c("forecasting_horizon",
                                   "model",
                                   "family"))


fwrite(out_selected,
       "C:\\Users\\gpitt\\Documents\\Postdoc\\Torino\\Mortality\\results\\ita_validation_testing_selection_l2.csv")



