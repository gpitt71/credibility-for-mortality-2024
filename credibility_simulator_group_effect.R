# IBMPopSim ----
rm(list=ls())
Sys.setenv(LANGUAGE = "en")
library(IBMPopSim)
library(dplyr)
library(StMoMo)
library(ggplot2)

N <- 100000


Q1 <- 0.8
N1 <- round(N*Q1)
ages_1 <- -65:0
tmp_df <- data.frame("birth"=rep(ages_1,each=N),
                     "death"=rep(rep(NA,N),each=length(ages_1)),
                     "risk_cls"= 3)

ages <- -65:-50
age_classes <- length(ages)

pop_df <- data.frame("birth"=rep(ages,each=N-N1),
                     "death"=rep(rep(NA,N-N1),each=age_classes),
                     "risk_cls"= sample(c(1,2), 
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

# d_k <- apply(LCforecastMale$rates, 2, function(x) stepfun(66:100, x))
# breaks <- 1:29
# death_male <- piecewise_xy(breaks,d_k)

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
newI.risk_cls =1;
else
newI.risk_cls= 2;
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
  set.seed(1996)
  sim_out <- popsim(model = model,
                  initial_population = pop_init,
                  events_bounds = c('death' = death_max, "entry"=params$lambda, "exit"=max(params$mu)),
                  parameters = params,
                  time = t,
                  age_max = 110,
                  multithreading = TRUE)
  
  
  }

ages.fit <- 50:80
years.fit <- 0:(t-1)

N_groups <- 3

# find deaths and exposures data 

D1 <- death_table(sim_out$population[sim_out$population$risk_cls==1,],
                  ages = c(ages.fit,last(ages.fit)), period = c(years.fit,last(years.fit)))

D2<- death_table(sim_out$population[sim_out$population$risk_cls==2,],
                 ages = c(ages.fit,last(ages.fit)), period = c(years.fit,last(years.fit)))

D3<- death_table(sim_out$population[sim_out$population$risk_cls==3,],
                 ages = c(ages.fit,last(ages.fit)), period = c(years.fit,last(years.fit)))

E1 <- exposure_table(sim_out$population[sim_out$population$risk_cls==1,],
                     ages = c(ages.fit,last(ages.fit)), period = c(years.fit,last(years.fit)))

E2 <- exposure_table(sim_out$population[sim_out$population$risk_cls==2,],
                     ages = c(ages.fit,last(ages.fit)), period = c(years.fit,last(years.fit)))
E3 <- exposure_table(sim_out$population[sim_out$population$risk_cls==3,],
                     ages = c(ages.fit,last(ages.fit)), period = c(years.fit,last(years.fit)))

# create a super-population 

# Dtot <- EWStMoMoMale$Dxt[as.character(ages.fit),as.character(1950:(1950+max(years.fit)))]
# Etot <- EWStMoMoMale$Ext[as.character(ages.fit),as.character(1950:(1950+max(years.fit)))]

# colnames(Etot) <- colnames(Dtot) <- years.fit[1:(length(years.fit))]


datahat <- structure(list(Dxt = D1, 
                          Ext = E1, 
                          ages = ages.fit,
                          years = years.fit,
                          type = 'central', 
                          series = 'female', label = 'total'), 
                     class = "StMoMoData")


mortality_model_lc <- lc(link="log")
mortality_model_apc <- apc(link="log")
mortality_model_rh <- rh(link="log")

model_option<-"rh"


test_statistic <- nll0<- nll1 <- aic_m0<- aic_m1<- bic_m0 <- bic_m1<- NULL


# for(model_option in c("lc","apc","rh")){
  
  assign("mortality_model",get(paste0("mortality_model_",model_option)))


  list_of_extra_exposures <- list()

for (i in 2:N_groups) {
    tmp_list <- list()
  
  tmp_list[['Dxt']] <- get(paste0("D",i))
  tmp_list[['Ext']] <- get(paste0("E",i))

  list_of_extra_exposures[[i-1]] <- tmp_list
  
  }

if(model_option=="rh"){
  
  mortality_model_fit <- fit(mortality_model,
                                                  data = datahat,
                                                  years.fit = years.fit,
                                                  ages.fit = ages.fit,
                                                  start.ax = start.ax,
                                                  start.bx = start.bx,
                                                  start.gc = start.gc,
                                                  start.kt = start.kt,
                                                  list_of_extra_exposures = list_of_extra_exposures)
  
  }else{

    mortality_model_fit <- fit(mortality_model, 
                               data = datahat,
                               years.fit = years.fit,
                               ages.fit = ages.fit,
                               list_of_extra_exposures = list_of_extra_exposures)
}

if(model_option == "apc"){start.gc = mortality_model_fit$gc}
if(model_option == "lc"){
  start.ax = mortality_model_fit$ax
  start.bx = mortality_model_fit$bx
  start.kt = mortality_model_fit$kt
}
  
  
forecasting_horizon<-1

year.predict <- max(years.fit)+forecasting_horizon

cv.arima.kt <- auto.arima(as.numeric(mortality_model_fit$kt), ic="bic")

if(model_option != "lc"){
  cv.arima.gc <- auto.arima(as.numeric(mortality_model_fit$gc), ic="bic")
  gc.order <- unname(arimaorder(cv.arima.gc))
}else{
  
  gc.order <- c(1,1,0)
}

mortality_model_forecast <- forecast(mortality_model_fit, 
                                     kt.method = "iarima",
                                     gc.order=gc.order,
                                     kt.order=unname(arimaorder(cv.arima.kt)),
                                     h = forecasting_horizon)



muxt_hat <- mortality_model_forecast$fitted


C1 <- apply(D1,1,sum,na.rm=T)/apply(E1*muxt_hat,1,sum,na.rm=T)

C2 <- apply(D2,1,sum,na.rm=T)/apply(E2*muxt_hat,1,sum,na.rm=T)

C3 <- apply(D3,1,sum,na.rm=T)/apply(E3*muxt_hat,1,sum,na.rm=T)

Fxt_1 <- D1/E1
Fxt_2 <- D2/E2
Fxt_3 <- D3/E3


tmpvar1 <- ((Fxt_1-muxt_hat)^2)/(muxt_hat^2)
tmpvar2 <- ((Fxt_2-muxt_hat)^2)/(muxt_hat^2)
tmpvar3 <- ((Fxt_3-muxt_hat)^2)/(muxt_hat^2)
tmpvar1[tmpvar1==1]<- NA
tmpvar2[tmpvar2==1]<- NA
tmpvar3[tmpvar3==1]<- NA

varthetax_1 <- apply(tmpvar1,1,mean, na.rm=T)
varthetax_2 <- apply(tmpvar2,1,mean, na.rm=T)
varthetax_3 <- apply(tmpvar3,1,mean, na.rm=T)


Z_1 <- 1/(1+varthetax_1*apply(E1*muxt_hat,1,sum,na.rm=T))

Z_2 <- 1/(1+varthetax_2*apply(E2*muxt_hat,1,sum,na.rm=T))

Z_3 <- 1/(1+varthetax_3*apply(E3*muxt_hat,1,sum,na.rm=T))

1-Z_1
1-Z_2
1-Z_3

dt_Zs <- data.frame(
  Zs = c(1-Z_1,1-Z_2),
  label = rep(c("Z_1","Z_2"),each=length(Z_1)),
  ages.code = rep(as.integer(names(Z_2)),2))

text_size= 28
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


# ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\",model_option,"_simulation_zs.pdf"),
#        width = 8,
#        height= 5)


chosen_age="70"
Cs_frame <- data.frame(y=D1[chosen_age,]/(E1*muxt_hat)[chosen_age,],
                       x=as.integer(colnames(D1))) %>%
  filter(y>0)

Cs_frame %>%
  ggplot(aes(x=x,
             y=y)) +
  geom_point() +
  geom_hline(yintercept=C1[chosen_age], linetype = "dotted", size= 1.2, color="#a71429") +
  geom_hline(yintercept=C1[chosen_age]+sqrt(varthetax_1[chosen_age]), size= 1.2, color="#4169E1", linetype = "dotted") +
  geom_hline(yintercept=C1[chosen_age]-sqrt(varthetax_1[chosen_age]), size= 1.2, color="#4169E1", linetype = "dotted") +
  theme_bw() +
  scale_x_continuous(breaks=c(10,20,30,40,50))+
  theme(text = element_text(size = text_size),
        legend.position="top")+
  ylab("")+
  xlab("")+
  labs(color="")

# ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\",model_option,"_simulation_cs_1.pdf"),
#        width = 8,
#        height= 5)

Cs_frame <- data.frame(y=D2[chosen_age,]/(E2*muxt_hat)[chosen_age,],
                       x=as.integer(colnames(D2))) %>%
  filter(y>0)

Cs_frame %>%
  ggplot(aes(x=x,
             y=y)) +
  geom_point() +
  geom_hline(yintercept=C2[chosen_age], linetype = "dotted", size= 1.2, color="#a71429") +
  geom_hline(yintercept=C2[chosen_age]+sqrt(varthetax_2[chosen_age]), size= 1.2, color="#4169E1", linetype = "dotted") +
  geom_hline(yintercept=C2[chosen_age]-sqrt(varthetax_2[chosen_age]), size= 1.2, color="#4169E1", linetype = "dotted") +
  theme_bw() +
  scale_x_continuous(breaks=c(10,20,30,40,50))+
  theme(text = element_text(size = text_size),
        legend.position="top")+
  ylab("")+
  xlab("")+
  labs(color="")

# ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\",model_option,"_simulation_cs_2.pdf"),
#        width = 8,
#        height= 5)


C1_matrix <- matrix(rep(C1, dim(D1)[2]),
                    nrow = dim(D1)[1],
                    byrow= FALSE)

C2_matrix <- matrix(rep(C2, dim(D1)[2]),
                    nrow = dim(D1)[1],
                    byrow= FALSE)

lik_h0 <- sum(mortality_model_fit$wxt*(D1*log(E1*muxt_hat)-E1*muxt_hat-lfactorial(D1)),na.rm = TRUE)+sum(mortality_model_fit$wxt*(D2*log(E2*muxt_hat)-E2*muxt_hat-lfactorial(D2)),na.rm = TRUE)
lik_h1 <- sum(mortality_model_fit$wxt*(D1*log(E1*muxt_hat*C1_matrix)-E1*muxt_hat*C1_matrix-lfactorial(D1)),na.rm = TRUE)+sum(mortality_model_fit$wxt*(D2*log(E2*muxt_hat*C2_matrix)-E2*muxt_hat*C2_matrix-lfactorial(D2)),na.rm = TRUE)
lrt <- -2*(lik_h0- lik_h1)


test_statistic <- c(test_statistic,lrt)

nll0 <- c(nll0,-lik_h0)
nll1 <- c(nll1,-lik_h1)

npar <- length(coef(mortality_model_fit$fittingModel))
nobs <- sum(D1!=0)


aic_m0 <- c(aic_m0,
            2*(npar)-2*lik_h0)
aic_m1 <- c(aic_m1,
            2*(npar+N_groups*length(ages.fit))-2*lik_h1)
bic_m0 <- c(bic_m0,
            (npar)*log(nobs)-2*lik_h0)
bic_m1 <- c(bic_m1,
            (npar+N_groups*length(ages.fit))*log(nobs)-2*lik_h1)
# }


out <- data.frame(model=c("lc","apc","rh"),
                  test_statistic,
                  nll0,
                  nll1,
                  aic_m0,
                  aic_m1,
                  bic_m0,
                  bic_m1)


# save(out,
#      file= "C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\statistics_simulation.txt")
# 

