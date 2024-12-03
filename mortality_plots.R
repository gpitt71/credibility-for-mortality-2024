rm(list = ls())
library(StMoMo)
library(HMDHFDplus)
library(tidyr)
library(dplyr)
library(demography)
library(ggplot2)
library(RColorBrewer)
# Full model ----

G1.mx <- hmd.mx(country = "USA",
                "gabriele.pittarello@uniroma1.it",
                "pm869fqW8nxozS8!")

G2.mx <- hmd.mx(country = "USA",
                "gabriele.pittarello@uniroma1.it",
                "pm869fqW8nxozS8!")

ages.fit <- 50:89
years.fit <- 1951:2003

series_code = 'male'
series_code_2 = 'female'

mortality_model_lc <- lc(link = "log")
mortality_model_apc <- apc(link = "log")
mortality_model_rh <- rh(link = "log")

N_groups = 2


out <- NULL 

forecasting_horizon=15

for (model_option in c("lc", "apc", "rh")) {
    
    mse_0_val <- mse_1_val <- mse_2_val <- NULL
    mae_0_val <- mae_1_val <- mae_2_val <- NULL
    
    mse_0 <- mse_1 <- mse_2 <- NULL
    mae_0 <- mae_1 <- mae_2 <- NULL
    
    assign("mortality_model", get(paste0("mortality_model_", model_option)))
    
    
    E1 <- G1.mx$pop[[series_code]][as.character(ages.fit), as.character(years.fit)]
    E2 <- G2.mx$pop[[series_code_2]][as.character(ages.fit), as.character(years.fit)]
    
    D1 <- E1  * G1.mx$rate[[series_code]][as.character(ages.fit), as.character(years.fit)]
    D2 <- E2  * G2.mx$rate[[series_code_2]][as.character(ages.fit), as.character(years.fit)]
    
    Ext <- E1 #+ E2
    Dxt <-  D1 #+ D2
    
    datahat <- structure(list(Dxt = Dxt, 
                              Ext = Ext, 
                              ages = ages.fit, 
                              years = years.fit, 
                              type = 'central', 
                              series = series_code, label = 'total'), 
                         class = "StMoMoData")
    
    
    list_of_extra_exposures <- list()
    
    list_of_extra_exposures[[1]] <- list(Dxt=D2,
                                         Ext=E2)
    
    mortality_model_fit <- fit(mortality_model,
                               data = datahat,
                               years.fit = years.fit,
                               ages.fit = ages.fit,
                               list_of_extra_exposures = list_of_extra_exposures)
    
    
    cv.arima.kt <- auto.arima(as.numeric(mortality_model_fit$kt), ic ="bic")
    
    if(model_option != "lc"){
      cv.arima.gc <- auto.arima(as.numeric(mortality_model_fit$gc), ic ="bic")
      gc.order <- unname(arimaorder(cv.arima.gc))
    }else{
      
      gc.order <- c(1,1,0)
    }
    
    mortality_model_forecast <- forecast(mortality_model_fit, 
                                         kt.method = "iarima",
                                         gc.order=gc.order,
                                         kt.order=unname(arimaorder(cv.arima.kt)),
                                         h = forecasting_horizon)
    
    
    muxt_actual_1 <- G1.mx$rate[[series_code]][as.character(ages.fit), as.character(max(years.fit+1):max(years.fit+forecasting_horizon))]
    muxt_actual_2 <- G2.mx$rate[[series_code_2]][as.character(ages.fit), as.character(max(years.fit+1):max(years.fit+forecasting_horizon))]
  
    
    dt_old <- data.frame(
      period.code = rep((min(years.fit):max(years.fit)),2),
      rate.code = rep(c("actual_male","actual_female"),each=max(years.fit)-min(years.fit)+1),
      muxt= c(G1.mx$rate[[series_code]]["65", as.character(min(years.fit):max(years.fit))],
              G2.mx$rate[[series_code_2]]["65", as.character(min(years.fit):max(years.fit))])
    )
    
    dt <- data.frame(
      
      age.code = rep("65",15),
      period.code = max(years.fit+1):max(years.fit+forecasting_horizon),
      actual_male = muxt_actual_1["65", ],
      actual_female = muxt_actual_2["65", ],
      s1 = mortality_model_forecast$rates["65",]
      
      
    )
    
    
    
    
    muxt_hat_predicted <- mortality_model_forecast$rates
    
    muxt_hat <- mortality_model_forecast$fitted
    
    C1 <- apply(D1, 1, sum) / apply(E1 * muxt_hat, 1, sum)
    C1_mx <- matrix(rep(C1,forecasting_horizon),byrow = F,ncol = forecasting_horizon)
    C2 <- apply(D2, 1, sum) / apply(E2 * muxt_hat, 1, sum)
    C2_mx <- matrix(rep(C2,forecasting_horizon),byrow = F,ncol = forecasting_horizon)
    
    
    Fxt_1 <- D1 / E1
    Fxt_2 <- D2 / E2
    
    varthetax_1 <- apply(((Fxt_1 - muxt_hat) ^ 2) / (muxt_hat ^ 2), 1, mean)
    varthetax_2 <- apply(((Fxt_2 - muxt_hat) ^ 2) / (muxt_hat ^ 2), 1, mean)
    
    Z_1 <- 1 / (1 + varthetax_1 * apply(E1 * muxt_hat, 1, sum))
    Z_1_mx <- matrix(rep(Z_1,forecasting_horizon),byrow = F,ncol = forecasting_horizon)
    
    Z_2 <- 1 / (1 + varthetax_2 * apply(E2 * muxt_hat, 1, sum))
    Z_2_mx <- matrix(rep(Z_2,forecasting_horizon),byrow = F,ncol = forecasting_horizon)
    
    muhat1_mx <- Z_1_mx*muxt_hat_predicted + (1-Z_1_mx)* C1_mx * muxt_hat_predicted
    muhat2_mx <- Z_2_mx*muxt_hat_predicted + (1-Z_2_mx)* C2_mx * muxt_hat_predicted
    
    dt[["s2_male"]] <- muhat1_mx["65",]
    dt[["s2_female"]] <- muhat2_mx["65",]
    
    G1data <- StMoMoData(G1.mx, series = series_code, type = "central")
    mortality_model_fit1 <- fit(mortality_model,
                                data = G1data,
                                years.fit = years.fit,
                                ages.fit = ages.fit)
    
    cv.arima.kt_1 <- auto.arima(as.numeric(mortality_model_fit1$kt), ic ="bic")
    
    G2data <- StMoMoData(G2.mx, series = series_code_2, type = "central")
    mortality_model_fit2 <- fit(mortality_model,
                                data = G2data,
                                years.fit = years.fit,
                                ages.fit = ages.fit)
    
    cv.arima.kt_2 <- auto.arima(as.numeric(mortality_model_fit2$kt), ic ="bic")
    
    if(model_option != "lc"){
      cv.arima.gc_1 <- auto.arima(as.numeric(mortality_model_fit1$gc), ic ="bic")
      gc.order_1 <- unname(arimaorder(cv.arima.gc_1))
      cv.arima.gc_2 <- auto.arima(as.numeric(mortality_model_fit2$gc), ic ="bic")
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
    
    
    mm_forecast_1_mx <- mortality_model_forecast1$rates[as.character(ages.fit),as.character(max(years.fit+1):max(years.fit+forecasting_horizon))]
    mm_forecast_2_mx <- mortality_model_forecast2$rates[as.character(ages.fit),as.character(max(years.fit+1):max(years.fit+forecasting_horizon))]
    
    
    dt[["s3_male"]] <- mm_forecast_1_mx["65",]
    dt[["s3_female"]] <- mm_forecast_2_mx["65",]
    
    
    text_size= 28
    
    
    dt_cs <- data.frame(
      Cis = c(C1,C2),
      label = rep(c("C1","C2"),each=length(C1)),
      ages.code = rep(as.integer(rownames(mortality_model_forecast1$rates)),N_groups))
    
    dt_cs %>%
      ggplot(aes(x=ages.code,
                 y=Cis)) +
      geom_point(aes(colour = label)) +
      theme_bw() +
      scale_x_continuous(breaks=c(50,60,70,80,90))+
      theme(text = element_text(size = text_size),
            legend.position="top")+
      ylab("")+
      xlab("")+
      labs(color="")+
      scale_color_manual(values = c("C1" = "#4169E1", 
                                    "C2" = "#a71429"),
                         labels = c(expression(C[x]^1), expression(C[x]^2))) 
    
    ggsave(paste0("C:\\Users\\gpitt\\Pictures\\mortality_project\\",model_option,"Cs.pdf"),
           width = 8,
           height= 5)
    
    dt_Zs <- data.frame(
      Zs = c(1-Z_1,1-Z_2),
      label = rep(c("Z_1","Z_2"),each=length(Z_1)),
      ages.code = rep(as.integer(rownames(mortality_model_forecast1$rates)),N_groups))
    
    dt_Zs %>%
      ggplot(aes(x=ages.code,
                 y=Zs)) +
      geom_point(aes(colour = label)) +
      theme_bw() +
      scale_x_continuous(breaks=c(50,60,70,80,90))+
      theme(text = element_text(size = text_size),
            legend.position="top")+
      ylab("")+
      xlab("")+
      labs(color="")+
      scale_color_manual(values = c("Z_1" = "#4169E1", 
                                    "Z_2" = "#a71429"),
                         labels = c(expression(1-Z[x]^1), expression(1-Z[x]^2))) 
    
    
    ggsave(paste0("C:\\Users\\gpitt\\Pictures\\mortality_project\\",model_option,"Zs.pdf"),
           width = 8,
           height= 5)
    
    dt %>%
      select(-s2_female,-s2_male,-s3_male,-s3_female)%>%
      pivot_longer(cols=c("actual_male","actual_female","s1"),
                   names_to=c("rate.code"),
                   values_to=c("muxt")) %>%
      full_join(dt_old, by=c("period.code",
                             "rate.code",
                             "muxt")) %>%
      ggplot(aes(x=period.code,
                 y=log(muxt))) +
      geom_line(aes(colour = rate.code),lwd=1.1) +
      geom_vline(xintercept =2004) +
      theme_bw() +
      theme(text = element_text(size = text_size),
            legend.position="none")+
      ylab("")+
      xlab("")+
      labs(color="")+
      scale_color_manual(values = c("actual_male" = "#4169E1", 
                                    "actual_female" = "#a71429",
                                    "s1" = "#454555"),
                         labels = c("Observed Female", "Observed Male", "Predicted (S1)")) 
    
    
    ggsave(paste0("C:\\Users\\gpitt\\Pictures\\mortality_project\\",model_option,"s1.pdf"),
           width = 8,
           height= 5)
    
    
    
    dt %>%
      select(-s1,-s3_male,-s3_female)%>%
      pivot_longer(cols=c("actual_male","actual_female","s2_male","s2_female"),
                   names_to=c("rate.code"),
                   values_to=c("muxt")) %>%
      full_join(dt_old, by=c("period.code",
                             "rate.code",
                             "muxt")) %>%
      ggplot(aes(x=period.code,
                 y=log(muxt))) +
      geom_line(aes(colour = rate.code),lwd=1.1) +
      geom_vline(xintercept =2004) +
      theme_bw() +
      theme(text = element_text(size = text_size),
            legend.position="none")+
      ylab("")+
      xlab("")+
      labs(color="")+
      scale_color_manual(values = c("actual_male" = "#4169E1",
                                    "actual_female" = "#a71429",
                                    "s2_male" = "#636b2f",
                                    "s2_female" = "#FF6A7A"),
                         labels = c("Observed Female", "Observed Male", "Predicted (S2, Female)", "Predicted (S2, Male)"))
    
    
    ggsave(paste0("C:\\Users\\gpitt\\Pictures\\mortality_project\\",model_option,"s2.pdf"),
           width = 8,
           height= 5)
    
    
    
    dt %>%
      select(-s1,-s2_male,-s2_female)%>%
      pivot_longer(cols=c("actual_male","actual_female","s3_male","s3_female"),
                   names_to=c("rate.code"),
                   values_to=c("muxt")) %>%
      full_join(dt_old, by=c("period.code",
                             "rate.code",
                             "muxt")) %>%
      ggplot(aes(x=period.code,
                 y=log(muxt))) +
      geom_line(aes(colour = rate.code),lwd=1.1) +
      geom_vline(xintercept =2004) +
      theme_bw() +
      theme(text = element_text(size = text_size),
            legend.position="none")+
      ylab("")+
      xlab("")+
      labs(color="")+
      scale_color_manual(values = c("actual_male" = "#4169E1",
                                    "actual_female" = "#a71429",
                                    "s3_male" = brewer.pal(n = 8, name = "RdYlBu")[7],
                                    "s3_female" = brewer.pal(n = 8, name = "RdYlBu")[3]),
                         labels = c("Observed Female", "Observed Male", "Predicted (S3, Female)", "Predicted (S3, Male)"))
    
    
    ggsave(paste0("C:\\Users\\gpitt\\Pictures\\mortality_project\\",model_option,"s3.pdf"),
           width = 8,
           height= 5)
    
    # brewer.pal(n = 8, name = "RdYlBu")[3]
      
    #display.brewer.pal(n = 8, name = 'RdYlBu')
    
    
    }

