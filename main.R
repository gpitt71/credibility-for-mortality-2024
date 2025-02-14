# Import functionalities

source(
  "C:/Users/pwt887/Documents/GitHub/credibility-for-mortality-2024/data_generators.R"
)
source(
  "C:/Users/pwt887/Documents/GitHub/credibility-for-mortality-2024/data_preprocessing.R"
)
source(
  "C:/Users/pwt887/Documents/GitHub/credibility-for-mortality-2024/utils_credibility.R"
)


# Inspection of Z and theta on one model fit ----

ages_fit <- 15:85
years_fit<- 118:140
N_groups=3


dt_0 <- data_generator_hmd_lt(
  seed_input = 1,
  exposure_sup = 100000,
  exposure_sub = c(5000, 500)
)


dt_0 %>%
  filter(Age== 55) %>%
  pivot_longer(cols=c("p","p1","p2"),
               names_to = "group",
               values_to = "probabilities") %>%
  mutate(
    group=recode(group,
           p = "Group 0",
           p1 = "Group 1",
           p2 = "Group 2")) %>%
  ggplot(aes(x=Period, y=probabilities)) +
  geom_line(aes(colour = group), size = 1.5)+
  theme_bw()+
  theme(text = element_text(size = 28),
        legend.position = "top") +
  scale_color_manual(values = c("Group 0" = "#a71429", "Group 1" = "#4169E1", "Group 2" = "#2E8B57"))+
  ylab(expression(q[xt]^i)) +
  xlab("") +
  labs(color = "")

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","mortality_rates_simulation.pdf"),
       width = 8,
       height= 5)

data_pp <- data_preprocessing_scenario_2(
  data = dt_0,
  ages_fit = ages_fit,
  years_fit=years_fit

)


mortality_models_fit <- fit_mortality_models(data_pp, years_fit-min(years_fit), ages_fit,separate_exposures =F)


total_model <- fit_and_predict_total_model(
  data = dt_0,
  data_pp =  data_pp,
  scenario = 2,
  N_groups = N_groups,
  mortality_models_fit = mortality_models_fit,
  years_fit = years_fit-min(years_fit),
  ages_fit = ages_fit,
  forecasting_horizon = 1,
  bias=min(years_fit)
)

credibility_model <- fit_and_predict_credibility_models(
  data = dt_0,
  data_pp =  data_pp,
  scenario = 2,
  N_groups = N_groups,
  mortality_models_fit = mortality_models_fit,
  years_fit = years_fit-min(years_fit),
  ages_fit = ages_fit,
  forecasting_horizon = 1,
  bias=min(years_fit)
)


separate_model <- fit_and_predict_separate_models(
  data = dt_0,
  data_pp =  data_pp,
  scenario = 2,
  N_groups = N_groups,
  mortality_models_fit = mortality_models_fit,
  years_fit = years_fit-min(years_fit),
  ages_fit = ages_fit,
  forecasting_horizon = 1,
  bias=min(years_fit)
)


mmodels <- c("lc", "apc","rh")
llks <- c(mortality_models_fit$lc$loglik,
              mortality_models_fit$apc$loglik,
              mortality_models_fit$rh$loglik)
names(llks) <-mmodels
min(llks)


variance_plotter(credibility_model,
                 predictor = "lc")

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_variances.pdf"),
       width = 8,
       height= 5)

plot_predicted_event_trends(mortality_models_fit,
                            data_pp,
                            model_option ="lc")


ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_event_trends.pdf"),
       width = 8,
       height= 5)


p <- weights_plotter (credibility_model, N_groups, ages_fit, predictor = "apc",ages_breaks = c(25,50,75,85))

p +
  scale_color_manual(
    values = c("Z_1" = "#a71429","Z_2" = "#4169E1", "Z_3" = "#006400"),
    labels = c(expression(Z[x]^1), expression(Z[x]^2), expression(Z[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_simulation_zs.pdf"),
       width = 8,
       height= 5)


all_thetas_plotter(credibility_model, N_groups, ages_fit, predictor = "lc") +
  scale_color_manual(
    values = c("C1" = "#a71429","C2" = "#4169E1", "C3" = "#006400"),
    labels = c(expression(hat(theta)[x]^0), expression(hat(theta)[x]^1), expression(hat(theta)[x]^2))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_allthetas_zs.pdf"),
       width = 8,
       height= 5)

credibility_model$actual_data$full_pp_data$list_of_extra_exposures[[1]]$Ext
tail(credibility_model$actual_data$full_pp_data$list_of_extra_exposures[[2]]$Ext)



credibility_model$apc$varthetax_1
credibility_model$lc$Z_1
credibility_model$lc$C1
credibility_model$actual_data$muxt_actual_1
credibility_model$lc$varthetax_2
credibility_model$lc$Z_2
credibility_model$lc$C2
credibility_model$actual_data$muxt_actual_2
credibility_model$lc$varthetax_3
credibility_model$lc$Z_3
credibility_model$lc$C3
credibility_model$actual_data$muxt_actual_3

## LC

variance_plotter(credibility_model,
                 predictor = "apc")

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","apc","_variances.pdf"),
       width = 8,
       height= 5)

exposure_plotter(data_pp = data_pp,
                 predictor = "lc")

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","exposures_by_age.pdf"),
       width = 8,
       height= 5)

p <- weights_plotter (credibility_model, N_groups, ages_fit, predictor = "apc",ages_breaks = c(5,10,17,50, 60, 70, 80, 90))

p +
  scale_color_manual(
    values = c("Z_1" = "#a71429","Z_2" = "#4169E1", "Z_3" = "#006400"),
    labels = c(expression(Z[x]^1), expression(Z[x]^2), expression(Z[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_simulation_zs.pdf"),
       width = 8,
       height= 5)


# Plot of the msep ----
## LC ----

reference_age = 65 
forecasting_horizon=5
cv.arima.kt <- auto.arima(as.numeric(mortality_models_fit$lc$kt), ic = "bic")
reference_group = 2


mortality_model_forecast <- forecast(
  mortality_models_fit$lc,
  kt.method = "iarima",
  # gc.order = gc.order,
  kt.order = unname(arimaorder(cv.arima.kt)),
  h = forecasting_horizon
)


ll <- dim(mortality_model_forecast$rates)[2]

muxt_hat_predicted <- mortality_model_forecast$rates[reference_age-min(ages_fit),(ll-forecasting_horizon):ll]

sigmabar_sim <- simulate(mortality_models_fit$lc, nsim = 1000, h = forecasting_horizon)

sigmabar <- apply(sigmabar_sim$rates[reference_age-min(ages_fit),,],1,var)

varthetax<- credibility_model$lc$varthetax_2[reference_age-min(ages_fit)]

term1 <- (sigmabar+muxt_hat_predicted^2)*(varthetax)+sigmabar

tmp_exp <- data_pp$list_of_extra_exposures[[reference_group-1]]$Ext[reference_age-min(ages_fit),]

fitted_muxt <- mortality_model_forecast$fitted[reference_age-min(ages_fit),]

term2 <- varthetax*(sum((tmp_exp^2) *( fitted_muxt)^2)/ sum(tmp_exp * fitted_muxt)^2)+(1/ sum(tmp_exp * fitted_muxt))

msep <- term1+term2

zed <- 1 / (1 + varthetax * sum(tmp_exp * fitted_muxt)) #1-z in papers notation...

msep <- term1+term2*(muxt_hat_predicted^2)*(1-zed)^2

thetaxi <- credibility_model$lc$C2[reference_age-min(ages_fit)]

predicted_rates <- zed*muxt_hat_predicted+(1-zed)*thetaxi*muxt_hat_predicted

data_pp2 <- data_preprocessing_scenario_2(
  data = dt_0,
  ages_fit = ages_fit,
  years_fit=c(years_fit,(last(years_fit)+1):(last(years_fit)+forecasting_horizon))
  
)

true_rates <- data_pp2$list_of_extra_exposures[[reference_group-1]]$Dxt[reference_age-min(ages_fit),]/data_pp2$list_of_extra_exposures[[reference_group-1]]$Ext[reference_age-min(ages_fit),]

zeros <- rep(NA,length(true_rates)-forecasting_horizon)
zeros_b <- rep(NA,length(true_rates)-forecasting_horizon-1)

sep_model_fit <- fit(
  lc(link="log"),
  data = structure(
    list(
      Dxt = data_pp$list_of_extra_exposures[[reference_group - 1]]$Dxt[as.character(ages_fit), as.character(years_fit-min(years_fit))],
      Ext = data_pp$list_of_extra_exposures[[reference_group - 1]]$Ext[as.character(ages_fit), as.character(years_fit-min(years_fit))],
      ages = ages_fit,
      years = years_fit,
      type = 'central',
      series = 'male',
      label =  paste0("superpop", reference_group)
    ),
    class = "StMoMoData"
  ),
  years.fit = years_fit,
  ages.fit = ages_fit,
  verbose = FALSE
)

cv.arima.kt <- auto.arima(as.numeric(sep_model_fit$kt), ic = "bic")

mortality_model_forecast_separate <- forecast(
  sep_model_fit,
  kt.method = "iarima",
  # gc.order = gc.order,
  kt.order = unname(arimaorder(cv.arima.kt)),
  h = forecasting_horizon
)

muxt_hat_predicted_separate <- mortality_model_forecast_separate$rates[reference_age-min(ages_fit),(ll-forecasting_horizon):ll]


# separate_bootstrap_model <- bootstrap(sep_model_fit, nBoot = 500, type = "residual")

# sep_boot_sim<- simulate(separate_bootstrap_model, h = 5)

# save.image("C:/Users/pwt887/Documents/GitHub/credibility-for-mortality-2024/output/my_environment.RData")

load("C:/Users/pwt887/Documents/GitHub/credibility-for-mortality-2024/output/my_environment.RData")

mxtvar <- apply(sep_boot_sim$rates, c(1, 2), var)

mxtvar <-mxtvar[reference_age-min(ages_fit),]

dt_all_rates <- data.frame(
  true_rates=true_rates,
  credibility_rates = c(zeros,predicted_rates),
  separate_rates = c(zeros,muxt_hat_predicted_separate),
  conf1_cred = c(zeros,predicted_rates+sqrt(msep)),
  conf2_cred = c(zeros,predicted_rates-sqrt(msep)),
  conf1_sep = c(zeros,muxt_hat_predicted_separate+sqrt(mxtvar)),
  conf2_sep = c(zeros,muxt_hat_predicted_separate-sqrt(mxtvar))  
)

dt_all_rates <- apply(dt_all_rates,2,log)

# library(tidyr)
# library(dplyr)

# Assume your dataframe is called dt_all_rates
dt_long <- dt_all_rates %>%  as.data.frame()%>%
  mutate(calendar_time=0:(length(true_rates)-1)) %>%
  as.data.frame()%>%
  pivot_longer(cols = c(true_rates, credibility_rates, separate_rates, 
                        conf1_cred, conf2_cred, conf1_sep, conf2_sep),
               names_to = "quantity",
               values_to = "values")

dt_long[['ltype']] <- "solid"
dt_long$ltype[grepl("conf",dt_long$quantity)] <- "dotted"


ggplot(dt_long, aes(x = calendar_time, y = values, color = quantity)) +
  geom_line(aes(linetype = ltype),
            size = 1) +
  scale_linetype_manual(values = c("solid" = "solid", "dotted" = "dotted")) +
  scale_color_manual(values = c("true_rates" = "black", 
                                "credibility_rates" = "#4169E1", 
                                "separate_rates" = "#a71429", 
                                "conf1_cred" = "#4169E1", 
                                "conf2_cred" = "#4169E1", 
                                "conf1_sep" = "#a71429", 
                                "conf2_sep" = "#a71429")) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "",
       color = "") +
  theme(legend.position = "none")

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","msep_65_g2.pdf"),
       width = 8,
       height= 5)

## Reference group 3 ----
reference_age=64
reference_group = 3


mortality_model_forecast <- forecast(
  mortality_models_fit$lc,
  kt.method = "iarima",
  # gc.order = gc.order,
  kt.order = unname(arimaorder(cv.arima.kt)),
  h = forecasting_horizon
)


ll <- dim(mortality_model_forecast$rates)[2]

muxt_hat_predicted <- mortality_model_forecast$rates[reference_age-min(ages_fit),(ll-forecasting_horizon):ll]

sigmabar_sim <- simulate(mortality_models_fit$lc, nsim = 1000, h = forecasting_horizon)

sigmabar <- apply(sigmabar_sim$rates[reference_age-min(ages_fit),,],1,var)

varthetax<- credibility_model$lc$varthetax_3[reference_age-min(ages_fit)]

term1 <- (sigmabar+muxt_hat_predicted^2)*(varthetax)+sigmabar

tmp_exp <- data_pp$list_of_extra_exposures[[reference_group-1]]$Ext[reference_age-min(ages_fit),]

fitted_muxt <- mortality_model_forecast$fitted[reference_age-min(ages_fit),]

term2 <- varthetax*(sum((tmp_exp^2) *( fitted_muxt)^2)/ sum(tmp_exp * fitted_muxt)^2)+(1/ sum(tmp_exp * fitted_muxt))

# msep <- term1+term2

zed <- 1 / (1 + varthetax * sum(tmp_exp * fitted_muxt)) #1-z in papers notation...

msep <- term1+term2*(muxt_hat_predicted^2)*(1-zed)^2

thetaxi <- credibility_model$lc$C2[reference_age-min(ages_fit)]

predicted_rates <- zed*muxt_hat_predicted+(1-zed)*thetaxi*muxt_hat_predicted

data_pp2 <- data_preprocessing_scenario_2(
  data = dt_0,
  ages_fit = ages_fit,
  years_fit=c(years_fit,(last(years_fit)+1):(last(years_fit)+forecasting_horizon))
  
)

true_rates <- data_pp2$list_of_extra_exposures[[reference_group-1]]$Dxt[reference_age-min(ages_fit),]/data_pp2$list_of_extra_exposures[[reference_group-1]]$Ext[reference_age-min(ages_fit),]

zeros <- rep(NA,length(true_rates)-forecasting_horizon)
zeros_b <- rep(NA,length(true_rates)-forecasting_horizon-1)

#algorithm fails we take global mttrend
sep_model_fit <- mortality_models_fit$lc

cv.arima.kt <- auto.arima(as.numeric(sep_model_fit$kt), ic = "bic")

mortality_model_forecast_separate <- forecast(
  sep_model_fit,
  kt.method = "iarima",
  # gc.order = gc.order,
  kt.order = unname(arimaorder(cv.arima.kt)),
  h = forecasting_horizon
)

muxt_hat_predicted_separate <- mortality_model_forecast_separate$rates[reference_age-min(ages_fit),(ll-forecasting_horizon):ll]


separate_bootstrap_model <- bootstrap(sep_model_fit, nBoot = 500, type = "residual")

sep_boot_sim<- simulate(separate_bootstrap_model, h = 5)

save.image("C:/Users/pwt887/Documents/GitHub/credibility-for-mortality-2024/output/my_environmentg3.RData")

load("C:/Users/pwt887/Documents/GitHub/credibility-for-mortality-2024/output/my_environmentg3.RData")

mxtvar <- apply(sep_boot_sim$rates, c(1, 2), var)

mxtvar <-mxtvar[reference_age-min(ages_fit),]

dt_all_rates <- data.frame(
  true_rates=true_rates,
  credibility_rates = c(zeros,predicted_rates),
  separate_rates = c(zeros,muxt_hat_predicted_separate),
  conf1_cred = c(zeros,predicted_rates+sqrt(msep)),
  conf2_cred = c(zeros,predicted_rates-sqrt(msep)),
  conf1_sep = c(zeros,muxt_hat_predicted_separate+sqrt(mxtvar)),
  conf2_sep = c(zeros,muxt_hat_predicted_separate-sqrt(mxtvar))  
)

dt_all_rates <- (apply(dt_all_rates,2,log))



# library(tidyr)
# library(dplyr)

# Assume your dataframe is called dt_all_rates
dt_long <- dt_all_rates %>%
  as.data.frame()%>%
  mutate(calendar_time=0:(length(true_rates)-1)) %>%
  pivot_longer(cols = c(true_rates, credibility_rates, separate_rates, 
                        conf1_cred, conf2_cred, conf1_sep, conf2_sep),
               names_to = "quantity",
               values_to = "values")

dt_long[['ltype']] <- "solid"
dt_long$ltype[grepl("conf",dt_long$quantity)] <- "dotted"


ggplot(dt_long, aes(x = calendar_time, y = values, color = quantity)) +
  geom_line(aes(linetype = ltype),
            size = 1) +
  scale_linetype_manual(values = c("solid" = "solid", "dotted" = "dotted")) +
  scale_color_manual(values = c("true_rates" = "black", 
                                "credibility_rates" = "#4169E1", 
                                "separate_rates" = "#a71429", 
                                "conf1_cred" = "#4169E1", 
                                "conf2_cred" = "#4169E1", 
                                "conf1_sep" = "#a71429", 
                                "conf2_sep" = "#a71429")) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "",
       color = "") +
  theme(legend.position = "none")

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","msep_65_g3.pdf"),
       width = 8,
       height= 5)

## Reference group 1 ----
reference_age=55
reference_group = 1


mortality_model_forecast <- forecast(
  mortality_models_fit$lc,
  kt.method = "iarima",
  # gc.order = gc.order,
  kt.order = unname(arimaorder(cv.arima.kt)),
  h = forecasting_horizon
)


ll <- dim(mortality_model_forecast$rates)[2]

muxt_hat_predicted <- mortality_model_forecast$rates[reference_age-min(ages_fit),(ll-forecasting_horizon):ll]

sigmabar_sim <- simulate(mortality_models_fit$lc, nsim = 1000, h = forecasting_horizon)

sigmabar <- apply(sigmabar_sim$rates[reference_age-min(ages_fit),,],1,var)

varthetax<- credibility_model$lc$varthetax_1[reference_age-min(ages_fit)]

term1 <- (sigmabar+muxt_hat_predicted^2)*(varthetax)+sigmabar

tmp_exp <- data_pp$datahat$Ext[reference_age-min(ages_fit),]

fitted_muxt <- mortality_model_forecast$fitted[reference_age-min(ages_fit),]

term2 <- varthetax*(sum((tmp_exp^2) *( fitted_muxt)^2)/ sum(tmp_exp * fitted_muxt)^2)+(1/ sum(tmp_exp * fitted_muxt))

# msep <- term1+term2

zed <- 1 / (1 + varthetax * sum(tmp_exp * fitted_muxt)) #1-z in papers notation...

msep <- term1+term2*(muxt_hat_predicted^2)*(1-zed)^2

thetaxi <- credibility_model$lc$C2[reference_age-min(ages_fit)]

predicted_rates <- zed*muxt_hat_predicted+(1-zed)*thetaxi*muxt_hat_predicted

data_pp2 <- data_preprocessing_scenario_2(
  data = dt_0,
  ages_fit = ages_fit,
  years_fit=c(years_fit,(last(years_fit)+1):(last(years_fit)+forecasting_horizon))
  
)

true_rates <- data_pp2$datahat$Dxt[reference_age-min(ages_fit),]/data_pp2$datahat$Ext[reference_age-min(ages_fit),]

zeros <- rep(NA,length(true_rates)-forecasting_horizon)
zeros_b <- rep(NA,length(true_rates)-forecasting_horizon-1)


sep_model_fit <- fit(
  lc(link="log"),
  data = structure(
    list(
      Dxt = data_pp$datahat$Dxt[as.character(ages_fit), as.character(years_fit-min(years_fit))],
      Ext = data_pp$datahat$Ext[as.character(ages_fit), as.character(years_fit-min(years_fit))],
      ages = ages_fit,
      years = years_fit,
      type = 'central',
      series = 'male',
      label =  paste0("superpop", reference_group)
    ),
    class = "StMoMoData"
  ),
  years.fit = years_fit,
  ages.fit = ages_fit,
  verbose = FALSE
)

cv.arima.kt <- auto.arima(as.numeric(sep_model_fit$kt), ic = "bic")

mortality_model_forecast_separate <- forecast(
  sep_model_fit,
  kt.method = "iarima",
  # gc.order = gc.order,
  kt.order = unname(arimaorder(cv.arima.kt)),
  h = forecasting_horizon
)

muxt_hat_predicted_separate <- mortality_model_forecast_separate$rates[reference_age-min(ages_fit),(ll-forecasting_horizon):ll]


separate_bootstrap_model <- bootstrap(sep_model_fit, nBoot = 500, type = "residual")

sep_boot_sim<- simulate(separate_bootstrap_model, h = 5)

save.image("C:/Users/pwt887/Documents/GitHub/credibility-for-mortality-2024/output/my_environmentg1.RData")


mxtvar <- apply(sep_boot_sim$rates, c(1, 2), var)

mxtvar <-mxtvar[reference_age-min(ages_fit),]

dt_all_rates <- data.frame(
  true_rates=true_rates,
  credibility_rates = c(zeros,predicted_rates),
  separate_rates = c(zeros,muxt_hat_predicted_separate),
  conf1_cred = c(zeros,predicted_rates+sqrt(msep)),
  conf2_cred = c(zeros,predicted_rates-sqrt(msep)),
  conf1_sep = c(zeros,muxt_hat_predicted_separate+sqrt(mxtvar)),
  conf2_sep = c(zeros,muxt_hat_predicted_separate-sqrt(mxtvar))  
)

dt_all_rates <- apply(dt_all_rates,2,log)

dt_all_rates[['calendar_time']] <- 0:(length(true_rates)-1) 

# library(tidyr)
# library(dplyr)

# Assume your dataframe is called dt_all_rates
dt_long <- dt_all_rates %>%
  as.data.frame()%>%
  pivot_longer(cols = c(true_rates, credibility_rates, separate_rates, 
                        conf1_cred, conf2_cred, conf1_sep, conf2_sep),
               names_to = "quantity",
               values_to = "values")

dt_long[['ltype']] <- "solid"
dt_long$ltype[grepl("conf",dt_long$quantity)] <- "dotted"


ggplot(dt_long, aes(x = calendar_time, y = values, color = quantity)) +
  geom_line(aes(linetype = ltype),
            size = 1) +
  scale_linetype_manual(values = c("solid" = "solid", "dotted" = "dotted")) +
  scale_color_manual(values = c("true_rates" = "black", 
                                "credibility_rates" = "#4169E1", 
                                "separate_rates" = "#a71429", 
                                "conf1_cred" = "#4169E1", 
                                "conf2_cred" = "#4169E1", 
                                "conf1_sep" = "#a71429", 
                                "conf2_sep" = "#a71429")) +
  theme_minimal() +
  labs(title = "",
       x = "",
       y = "",
       color = "") +
  theme(legend.position = "none")

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","msep_65_g1.pdf"),
       width = 8,
       height= 5)

# Assess the model performance in age breakets ----

ages_fit <- 15:85
years_fit_basic <- 118:140
N_groups=3
forecasting_horizon=1

ageb_width <- 10
out <- NULL

for (seed_ix in 1:3) {

  dt_0 <- data_generator_hmd_lt(
    seed_input = seed_ix,
    exposure_sup = 100000,
    exposure_sub = c(5000, 500)
  )

    for (years_ix in 1:6) {


      years_fit <- c(years_fit_basic, (last(years_fit_basic)+1):(last(years_fit_basic)+years_ix))

      rm(data_pp,mortality_models_fit,total_model,credibility_model,separate_model,performance_total_df,performance_credibility_df,performance_separate_df)

      data_pp <- data_preprocessing_scenario_2(
        data = dt_0,
        ages_fit = ages_fit,
        years_fit=years_fit

      )

      mortality_models_fit <- fit_mortality_models(data_pp, years_fit-min(years_fit), ages_fit,separate_exposures =F)


      total_model <- fit_and_predict_total_model(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 2,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit-min(years_fit),
        ages_fit = ages_fit,
        forecasting_horizon = 1,
        bias=min(years_fit)
      )

      credibility_model <- fit_and_predict_credibility_models(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 2,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit-min(years_fit),
        ages_fit = ages_fit,
        forecasting_horizon = 1,
        bias=min(years_fit)
      )


      separate_model <- fit_and_predict_separate_models(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 2,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit-min(years_fit),
        ages_fit = ages_fit,
        forecasting_horizon = 1,
        bias=min(years_fit)
      )

      tot_ages <- max(ages_fit)-min(ages_fit)

      tot_intervals <- tot_ages/ageb_width
      answer_total=NULL
      answer_credibility=NULL
      answer_mle=NULL
      answer_separate=NULL

      for(pfix in 1:tot_intervals){

        subset_rates <- (min(ages_fit)+(pfix-1)*ageb_width):(min(ages_fit)+(pfix)*ageb_width)
        interval_name <- paste(min(ages_fit)+(pfix-1)*ageb_width+1,min(ages_fit)+(pfix)*ageb_width,sep = "-")

        performance_total <- model_assessment(
          model_fit_and_prediction = total_model,
          N_groups = N_groups,
          years_fit = years_fit,
          forecast = forecasting_horizon,
          subset_rates = subset_rates,
          bias=min(ages_fit)
        )

        performance_total_df <- do.call(rbind,
                                        lapply(performance_total, data.frame, stringsAsFactors = FALSE))
        performance_total_df[['seed']] <- seed_ix
        performance_total_df[['model']] <- 'total'
        performance_total_df[['predictor']] <- rownames(performance_total_df)
        performance_total_df[['years_ix']] <- years_ix
        performance_total_df[['forecasting_horizon']] <- forecasting_horizon
        performance_total_df[['age_breakets']] <- interval_name
        rownames(performance_total_df)<-NULL

        answer_total<-rbind(answer_total,performance_total_df)


        performance_credibility <- model_assessment(
          model_fit_and_prediction = credibility_model,
          N_groups = N_groups,
          years_fit = years_fit,
          forecast = forecasting_horizon,
          subset_rates = subset_rates,
          bias=min(ages_fit)
        )

        performance_credibility_df <- do.call(rbind,
                                              lapply(performance_credibility, data.frame, stringsAsFactors = FALSE))
        performance_credibility_df[['seed']] <- seed_ix
        performance_credibility_df[['model']] <- 'credibility'
        performance_credibility_df[['years_ix']] <- years_ix
        performance_credibility_df[['forecasting_horizon']] <- forecasting_horizon
        performance_credibility_df[['predictor']] <- rownames(performance_credibility_df)
        performance_credibility_df[['age_breakets']] <- interval_name
        rownames(performance_credibility_df)<-NULL

        answer_credibility<-rbind(answer_credibility,performance_credibility_df)

        performance_credibility_full_mle <-model_assessment_full_mle(credibility_model,
                                                                     N_groups,
                                                                     years_fit,
                                                                     forecasting_horizon,
                                                                     subset_rates = subset_rates,
                                                                     bias=min(ages_fit))


        performance_credibility_full_mle_df <- do.call(rbind,
                                                       lapply(performance_credibility_full_mle, data.frame, stringsAsFactors = FALSE))
        performance_credibility_full_mle_df[['seed']] <- seed_ix
        performance_credibility_full_mle_df[['model']] <- 'full_mle'
        performance_credibility_full_mle_df[['years_ix']] <- years_ix
        performance_credibility_full_mle_df[['forecasting_horizon']] <- forecasting_horizon
        performance_credibility_full_mle_df[['predictor']] <- rownames(performance_credibility_full_mle_df)
        performance_credibility_full_mle_df[['age_breakets']] <- interval_name
        rownames(performance_credibility_full_mle_df)<-NULL

        answer_mle<-rbind(answer_mle,performance_credibility_full_mle_df)

        performance_separate <- model_assessment(
          model_fit_and_prediction = separate_model,
          N_groups = N_groups,
          years_fit = years_fit,
          forecast = forecasting_horizon,
          subset_rates = subset_rates,
          bias=min(ages_fit)
        )

        performance_separate_df <- do.call(rbind,
                                           lapply(performance_separate, data.frame, stringsAsFactors = FALSE))
        performance_separate_df[['seed']] <- seed_ix
        performance_separate_df[['model']] <- 'separate'
        performance_separate_df[['years_ix']] <- years_ix
        performance_separate_df[['forecasting_horizon']] <- forecasting_horizon
        performance_separate_df[['predictor']] <- rownames(performance_separate_df)
        performance_separate_df[['age_breakets']] <- interval_name
        rownames(performance_separate_df)<-NULL

        answer_separate<-rbind(answer_separate,performance_separate_df)

      }


      out <- rbind(
        out,
        answer_total,
        answer_credibility,
        answer_mle,
        answer_separate
      )


    }

}


data.table::fwrite(out,"C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\results_age_breakets.txt")





## LC: thetas 75----

# years_fit = 0:29

thetas_plotter(
  credibility_model,
  subgroup=1,
  ages_fit,
  years_fit,
  mortality_models_fit,
  predictor = "lc",
  chosen_age = 60
)

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_simulation_cs_1.pdf"),
       width = 8,
       height= 5)

thetas_plotter(
  credibility_model,
  subgroup=2,
  ages_fit,
  years_fit,
  mortality_models_fit,
  predictor = "lc",
  chosen_age = 60
)

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_simulation_cs_2.pdf"),
       width = 8,
       height= 5)

thetas_plotter(
  credibility_model,
  subgroup=3,
  ages_fit,
  years_fit,
  mortality_models_fit,
  predictor = "lc",
  chosen_age = 59
)

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_simulation_cs_3.pdf"),
       width = 8,
       height= 5)

## LC: thetas ----

all_thetas_plotter(credibility_model, N_groups, ages_fit, predictor = "lc") +
  scale_color_manual(
    values = c("C1" = "#a71429","C2" = "#4169E1", "C3" = "#006400"),
    labels = c(expression(hat(theta)[x]^1), expression(hat(theta)[x]^2), expression(hat(theta)[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_allthetas_zs.pdf"),
       width = 8,
       height= 5)

## APC: thetas ----

all_thetas_plotter(credibility_model, N_groups, ages_fit, predictor = "apc") +
  scale_color_manual(
    values = c("C1" = "#a71429","C2" = "#4169E1", "C3" = "#006400"),
    labels = c(expression(hat(Theta)[x]^1), expression(hat(Theta)[x]^2), expression(hat(Theta)[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","apc","_allthetas_zs.pdf"),
       width = 8,
       height= 5)

## RH: thetas ----

all_thetas_plotter(credibility_model, N_groups, ages_fit, predictor = "rh") +
  scale_color_manual(
    values = c("C1" = "#a71429","C2" = "#4169E1", "C3" = "#006400"),
    labels = c(expression(hat(Theta)[x]^1), expression(hat(Theta)[x]^2), expression(hat(Theta)[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","rh","_allthetas_zs.pdf"),
       width = 8,
       height= 5)

## LC: weights ----

p <- weights_plotter (credibility_model, N_groups, ages_fit, predictor = "lc")

p +
  scale_color_manual(
    values = c("Z_1" = "#a71429","Z_2" = "#4169E1", "Z_3" = "#006400"),
    labels = c(expression(Z[x]^1), expression(Z[x]^2), expression(Z[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_simulation_zs.pdf"),
       width = 8,
       height= 5)

## APC: weights ----

p <- weights_plotter (credibility_model, N_groups, ages_fit, predictor = "apc")

p +
  scale_color_manual(
    values = c("Z_1" = "#a71429","Z_2" = "#4169E1", "Z_3" = "#006400"),
    labels = c(expression(Z[x]^1), expression(Z[x]^2), expression(Z[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","apc","_simulation_zs.pdf"),
       width = 8,
       height= 5)

## RH: weights ----

p <- weights_plotter (credibility_model, N_groups, ages_fit, predictor = "rh")

p +
  scale_color_manual(
    values = c("Z_1" = "#a71429","Z_2" = "#4169E1", "Z_3" = "#006400"),
    labels = c(expression(Z[x]^1), expression(Z[x]^2), expression(Z[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","rh","_simulation_zs.pdf"),
       width = 8,
       height= 5)





# Some parameters ----
ages_fit <- 50:65
years_fit_basic <- 0:28
#
# years_fit <- 0:29
N_groups <- 3
# forecasting_horizon <-1
# years_ix<-1


# Simulate data: scenario 1 ----

out <- NULL

for (seed_ix in 1:3) {#1:3
  dt_0 <- data_generator_1(seed = seed_ix)

  for (forecasting_horizon in c(1,5,8)) {
    # 50 years is the maximum amount of data that we have

    for (years_ix in 1:6) {
      years_fit <- c(years_fit_basic, (last(years_fit_basic)+1):(last(years_fit_basic)+years_ix))

      rm(data_pp,mortality_models_fit,total_model,credibility_model,separate_model,performance_total_df,performance_credibility_df,performance_separate_df)

      data_pp <- data_preprocessing_scenario_1(
        data = dt_0,
        N_groups = N_groups,
        ages_fit = ages_fit,
        years_fit = years_fit
      )

      mortality_models_fit <- fit_mortality_models(data_pp, years_fit, ages_fit)

      total_model <- fit_and_predict_total_model(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 1,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit,
        ages_fit = ages_fit,
        forecasting_horizon = forecasting_horizon
      )

      credibility_model <- fit_and_predict_credibility_models(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 1,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit,
        ages_fit = ages_fit,
        forecasting_horizon = forecasting_horizon
      )



      separate_model <- fit_and_predict_separate_models(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 1,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit,
        ages_fit = ages_fit,
        forecasting_horizon = forecasting_horizon
      )


      performance_total <- model_assessment(
        model_fit_and_prediction = total_model,
        N_groups = N_groups,
        years_fit = years_fit,
        forecast = forecasting_horizon
      )

      performance_total_df <- do.call(rbind,
                                      lapply(performance_total, data.frame, stringsAsFactors = FALSE))
      performance_total_df[['seed']] <- seed_ix
      performance_total_df[['model']] <- 'total'
      performance_total_df[['predictor']] <- rownames(performance_total_df)
      performance_total_df[['years_ix']] <- years_ix
      performance_total_df[['forecasting_horizon']] <- forecasting_horizon
      rownames(performance_total_df)<-NULL

      performance_credibility <- model_assessment(
        model_fit_and_prediction = credibility_model,
        N_groups = N_groups,
        years_fit = years_fit,
        forecast = forecasting_horizon
      )

      performance_credibility_df <- do.call(rbind,
                                            lapply(performance_credibility, data.frame, stringsAsFactors = FALSE))
      performance_credibility_df[['seed']] <- seed_ix
      performance_credibility_df[['model']] <- 'credibility'
      performance_credibility_df[['years_ix']] <- years_ix
      performance_credibility_df[['forecasting_horizon']] <- forecasting_horizon
      performance_credibility_df[['predictor']] <- rownames(performance_credibility_df)
      rownames(performance_credibility_df)<-NULL


      performance_separate <- model_assessment(
        model_fit_and_prediction = separate_model,
        N_groups = N_groups,
        years_fit = years_fit,
        forecast = forecasting_horizon
      )

      performance_separate_df <- do.call(rbind,
                                         lapply(performance_separate, data.frame, stringsAsFactors = FALSE))
      performance_separate_df[['seed']] <- seed_ix
      performance_separate_df[['model']] <- 'separate'
      performance_separate_df[['years_ix']] <- years_ix
      performance_separate_df[['forecasting_horizon']] <- forecasting_horizon
      performance_separate_df[['predictor']] <- rownames(performance_separate_df)
      rownames(performance_separate_df)<-NULL

      out <- rbind(
        out,
        performance_total_df,
        performance_credibility_df,
        performance_separate_df
      )

    }
  }

}

out[['scenario']] <- 1

data.table::fwrite(out,"C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\results_scenario1.txt")





# Simulate data: scenario 2 ----
N_groups=3
out <- NULL
ages_fit <- 50:65
years_fit_basic <- 0:15

for (seed_ix in 1:3) {#1:3

  dt_0 <- data_generator_hmd(seed_input=seed_ix,
                             exposure_sup=100000,
                             exposure_sub=c(5000,
                                            500))

  for (forecasting_horizon in c(1,5,8)) {
    # 50 years is the maximum amount of data that we have

    for (years_ix in 1:6) {
      years_fit <- c(years_fit_basic, (last(years_fit_basic)+1):(last(years_fit_basic)+years_ix))

      rm(data_pp,mortality_models_fit,total_model,credibility_model,separate_model,performance_total_df,performance_credibility_df,performance_separate_df)


      data_pp <- data_preprocessing_scenario_2(
        data = dt_0,
        ages_fit = ages_fit,
        years_fit=70:(70+max(years_fit))

      )

      # data_pp <- data_preprocessing_scenario_1(
      #   data = dt_0,
      #   N_groups = N_groups,
      #   ages_fit = ages_fit,
      #   years_fit = years_fit
      # )

      mortality_models_fit <- fit_mortality_models(data_pp,
                                                   years_fit,
                                                   ages_fit,
                                                   separate_exposures=F)

      total_model <- fit_and_predict_total_model(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 2,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit,
        ages_fit = ages_fit,
        forecasting_horizon = forecasting_horizon,
        bias=min(70:100)
      )

      credibility_model <- fit_and_predict_credibility_models(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 2,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit,
        ages_fit = ages_fit,
        forecasting_horizon = forecasting_horizon,
        bias=min(70:100)
      )



      separate_model <- fit_and_predict_separate_models(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 2,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit,
        ages_fit = ages_fit,
        forecasting_horizon = forecasting_horizon,
        bias=min(70:100)
      )


      performance_total <- model_assessment(
        model_fit_and_prediction = total_model,
        N_groups = N_groups,
        years_fit = years_fit,
        forecast = forecasting_horizon
      )

      performance_total_df <- do.call(rbind,
                                      lapply(performance_total, data.frame, stringsAsFactors = FALSE))
      performance_total_df[['seed']] <- seed_ix
      performance_total_df[['model']] <- 'total'
      performance_total_df[['predictor']] <- rownames(performance_total_df)
      performance_total_df[['years_ix']] <- years_ix
      performance_total_df[['forecasting_horizon']] <- forecasting_horizon
      rownames(performance_total_df)<-NULL

      performance_credibility <- model_assessment(
        model_fit_and_prediction = credibility_model,
        N_groups = N_groups,
        years_fit = years_fit,
        forecast = forecasting_horizon
      )

      performance_credibility_df <- do.call(rbind,
                                            lapply(performance_credibility, data.frame, stringsAsFactors = FALSE))
      performance_credibility_df[['seed']] <- seed_ix
      performance_credibility_df[['model']] <- 'credibility'
      performance_credibility_df[['years_ix']] <- years_ix
      performance_credibility_df[['forecasting_horizon']] <- forecasting_horizon
      performance_credibility_df[['predictor']] <- rownames(performance_credibility_df)
      rownames(performance_credibility_df)<-NULL


      performance_credibility_full_mle <-model_assessment_full_mle(credibility_model,
                                N_groups,
                                years_fit,
                                forecasting_horizon)


      performance_credibility_full_mle_df <- do.call(rbind,
                                            lapply(performance_credibility_full_mle, data.frame, stringsAsFactors = FALSE))
      performance_credibility_full_mle_df[['seed']] <- seed_ix
      performance_credibility_full_mle_df[['model']] <- 'full_mle'
      performance_credibility_full_mle_df[['years_ix']] <- years_ix
      performance_credibility_full_mle_df[['forecasting_horizon']] <- forecasting_horizon
      performance_credibility_full_mle_df[['predictor']] <- rownames(performance_credibility_full_mle_df)
      rownames(performance_credibility_full_mle_df)<-NULL

      performance_separate <- model_assessment(
        model_fit_and_prediction = separate_model,
        N_groups = N_groups,
        years_fit = years_fit,
        forecast = forecasting_horizon
      )

      performance_separate_df <- do.call(rbind,
                                         lapply(performance_separate, data.frame, stringsAsFactors = FALSE))
      performance_separate_df[['seed']] <- seed_ix
      performance_separate_df[['model']] <- 'separate'
      performance_separate_df[['years_ix']] <- years_ix
      performance_separate_df[['forecasting_horizon']] <- forecasting_horizon
      performance_separate_df[['predictor']] <- rownames(performance_separate_df)
      rownames(performance_separate_df)<-NULL

      out <- rbind(
        out,
        performance_total_df,
        performance_credibility_df,
        performance_credibility_full_mle_df,
        performance_separate_df
      )

    }
  }

}

out[['scenario']] <- 2

data.table::fwrite(out,"C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\results_scenario2.txt")


# Simulate data: scenario 3 ----

out <- NULL
ages_fit <- 50:65
years_fit_basic <- 0:15

for (seed_ix in 1:3) {#1:3

  dt_0 <- data_generator_hmd_lt(
    seed_input = 1,
    exposure_sup = 100000,
    exposure_sub = c(5000, 500)
  )

  for (forecasting_horizon in c(1,5,8)) {
    # 50 years is the maximum amount of data that we have

    for (years_ix in 1:6) {
      years_fit <- c(years_fit_basic, (last(years_fit_basic)+1):(last(years_fit_basic)+years_ix))

      rm(data_pp,mortality_models_fit,total_model,credibility_model,separate_model,performance_total_df,performance_credibility_df,performance_separate_df)


      data_pp <- data_preprocessing_scenario_2(
        data = dt_0,
        ages_fit = ages_fit,
        years_fit=70:(70+max(years_fit))

      )

      # data_pp <- data_preprocessing_scenario_1(
      #   data = dt_0,
      #   N_groups = N_groups,
      #   ages_fit = ages_fit,
      #   years_fit = years_fit
      # )

      mortality_models_fit <- fit_mortality_models(data_pp,
                                                   years_fit,
                                                   ages_fit,
                                                   separate_exposures=T)

      total_model <- fit_and_predict_total_model(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 2,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit,
        ages_fit = ages_fit,
        forecasting_horizon = forecasting_horizon,
        bias=min(70:100)
      )

      credibility_model <- fit_and_predict_credibility_models(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 2,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit,
        ages_fit = ages_fit,
        forecasting_horizon = forecasting_horizon,
        bias=min(70:100)
      )



      separate_model <- fit_and_predict_separate_models(
        data = dt_0,
        data_pp =  data_pp,
        scenario = 2,
        N_groups = N_groups,
        mortality_models_fit = mortality_models_fit,
        years_fit = years_fit,
        ages_fit = ages_fit,
        forecasting_horizon = forecasting_horizon,
        bias=min(70:100)
      )


      performance_total <- model_assessment(
        model_fit_and_prediction = total_model,
        N_groups = N_groups,
        years_fit = years_fit,
        forecast = forecasting_horizon
      )

      performance_total_df <- do.call(rbind,
                                      lapply(performance_total, data.frame, stringsAsFactors = FALSE))
      performance_total_df[['seed']] <- seed_ix
      performance_total_df[['model']] <- 'total'
      performance_total_df[['predictor']] <- rownames(performance_total_df)
      performance_total_df[['years_ix']] <- years_ix
      performance_total_df[['forecasting_horizon']] <- forecasting_horizon
      rownames(performance_total_df)<-NULL

      performance_credibility <- model_assessment(
        model_fit_and_prediction = credibility_model,
        N_groups = N_groups,
        years_fit = years_fit,
        forecast = forecasting_horizon
      )

      performance_credibility_df <- do.call(rbind,
                                            lapply(performance_credibility, data.frame, stringsAsFactors = FALSE))
      performance_credibility_df[['seed']] <- seed_ix
      performance_credibility_df[['model']] <- 'credibility'
      performance_credibility_df[['years_ix']] <- years_ix
      performance_credibility_df[['forecasting_horizon']] <- forecasting_horizon
      performance_credibility_df[['predictor']] <- rownames(performance_credibility_df)
      rownames(performance_credibility_df)<-NULL


      performance_credibility_full_mle <-model_assessment_full_mle(credibility_model,
                                                                   N_groups,
                                                                   years_fit,
                                                                   forecasting_horizon)


      performance_credibility_full_mle_df <- do.call(rbind,
                                                     lapply(performance_credibility_full_mle, data.frame, stringsAsFactors = FALSE))
      performance_credibility_full_mle_df[['seed']] <- seed_ix
      performance_credibility_full_mle_df[['model']] <- 'full_mle'
      performance_credibility_full_mle_df[['years_ix']] <- years_ix
      performance_credibility_full_mle_df[['forecasting_horizon']] <- forecasting_horizon
      performance_credibility_full_mle_df[['predictor']] <- rownames(performance_credibility_full_mle_df)
      rownames(performance_credibility_full_mle_df)<-NULL

      performance_separate <- model_assessment(
        model_fit_and_prediction = separate_model,
        N_groups = N_groups,
        years_fit = years_fit,
        forecast = forecasting_horizon
      )

      performance_separate_df <- do.call(rbind,
                                         lapply(performance_separate, data.frame, stringsAsFactors = FALSE))
      performance_separate_df[['seed']] <- seed_ix
      performance_separate_df[['model']] <- 'separate'
      performance_separate_df[['years_ix']] <- years_ix
      performance_separate_df[['forecasting_horizon']] <- forecasting_horizon
      performance_separate_df[['predictor']] <- rownames(performance_separate_df)
      rownames(performance_separate_df)<-NULL

      out <- rbind(
        out,
        performance_total_df,
        performance_credibility_df,
        performance_credibility_full_mle_df,
        performance_separate_df
      )

    }
  }

}

out[['scenario']] <- 3

data.table::fwrite(out,"C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\results_scenario3.txt")










