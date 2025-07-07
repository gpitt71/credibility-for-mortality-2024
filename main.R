# Import functionalities

source(
  "~/credibility-for-mortality-2024/data_generators.R"
)

source(
  "~/credibility-for-mortality-2024/data_preprocessing.R"
)

source(
  "~/credibility-for-mortality-2024/utils_credibility.R"
)


# How to fit the different models ----

## Ages and years of interest
ages_fit <- 15:85
years_fit<- 118:140
N_groups=3

## Generate the data
dt_0 <- data_generator_hmd_lt(
  seed_input = 1,
  hmd_username="",
  hmd_password="",
  exposure_sup = 94500,
  exposure_sub = c(5000, 500)
)


## Preprocess the data
data_pp <- data_preprocessing_scenario_2(
  data = dt_0,
  ages_fit = ages_fit,
  years_fit=years_fit

)


## Global mortality trend fit
mortality_models_fit <- fit_mortality_models(data_pp, years_fit-min(years_fit), ages_fit,separate_exposures =F)



## Total mortality model - model D
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


## Credibility mortality model - model A (and B)
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


## Separate mortality model - model C
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


## Figure 6 ----

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
  geom_line(aes(colour = group), linewidth = 1.5)+
  theme_bw()+
  theme(text = element_text(size = 28),
        legend.position = "none") +
  scale_color_manual(values = c("Group 0" = "#a71429", "Group 1" = "#4169E1", "Group 2" = "#2E8B57"))+
  ylab(expression(q[xt]^i)) +
  xlab("") +
  labs(color = "")

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","mortality_rates_simulation.pdf"),
       width = 8,
       height= 5)

# Figure 2 ----
variance_plotter(credibility_model,
                 ages_fit=ages_fit,
                 predictor = "lc") + scale_fill_discrete(guide="none")

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_variances.pdf"),
       width = 8,
       height= 5)

plot_predicted_event_trends(mortality_models_fit,
                            data_pp,
                            model_option ="lc")

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_event_trends.pdf"),
       width = 8,
       height= 5)


p <- weights_plotter (credibility_model, N_groups, ages_fit, predictor = "lc",ages_breaks = c(25,50,75,85))


p +
  scale_color_manual(
    values = c("Z_2" = "#4169E1", "Z_3" = "#006400","Z_1" = "#a71429"),
    labels = c(expression(Z[x]^1), expression(Z[x]^2),expression(Z[x]^3))
  )

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_simulation_zs.pdf"),
       width = 8,
       height= 5)

#extract simulated age-group effects
dt<-print_age_effects(
  seed_input = 1,
  hmd_username="",
  hmd_password="",
  exposure_sup = 94500,
  exposure_sub = c(5000, 500)
) %>% as.data.frame()

dt[['ages']] <- 1:nrow(dt)

dt<-dt[ages_fit,]

all_thetas_plotter(credibility_model, N_groups, ages_fit, predictor = "lc") +
  scale_color_manual(
    values = c("C2" = "#4169E1", "C3" = "#006400","C1" = "#a71429"),
    labels = c( expression(hat(theta)[x]^1), expression(hat(theta)[x]^2),expression(hat(theta)[x]^3))
  )+ geom_line(data= dt, aes(x=ages,y=age_eff_2),size=1.5, alpha=.4, color="#006400")+
  geom_hline(yintercept=1.,size=1.5, alpha=.4, color="#a71429")+
  geom_line(data= dt, aes(x=ages,y=age_eff_1),size=1.5, alpha=.4, color="#4169E1")


ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_allthetas_zs.pdf"),
       width = 8,
       height= 5)

## The credibility model output can be seen here
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


# Section 4.4 ---- 
# Inspect the Mean Squared Error of Prediction (MSEP) ----
## Reference group 2 ----

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

zed <- 1 / (1 + varthetax * sum(tmp_exp * fitted_muxt)) #1-z in papers notation...

msep <- term1+term2*(muxt_hat_predicted^2)*(1-zed)^2

thetaxi <- credibility_model$lc$C2[reference_age-min(ages_fit)]

predicted_rates <- zed*muxt_hat_predicted+(1-zed)*thetaxi*muxt_hat_predicted

insample_rates <- zed*fitted_muxt+(1-zed)*thetaxi*fitted_muxt

our_rates <- c(insample_rates,
               predicted_rates)

tmp_exp_future <- data_pp$list_of_extra_exposures[[reference_group-1]]$Ext[reference_age-min(ages_fit),(ll-forecasting_horizon):ll]

{
  set.seed(1)
  simulation_of_observed_bundles <- sapply(c(tmp_exp,tmp_exp_future)*our_rates, function(x) rpois(1000, x))

  tmp_exp_matrix <- matrix(rep(c(tmp_exp,tmp_exp_future),
                               1000),
                           nrow = 1000,
                           byrow = T)


  simulation_of_observed_bundles <- simulation_of_observed_bundles/tmp_exp_matrix

  simulated_variance_if_real <- apply(simulation_of_observed_bundles, MARGIN = 2, sd)


  }




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

# save(sep_boot_sim, file = "~/credibility-for-mortality-2024/output/bootg2.RData")

load("~/credibility-for-mortality-2024/output/bootg2.RData")


# save.image("~/credibility-for-mortality-2024/output/my_environment.RData")

# load("~/credibility-for-mortality-2024/output/my_environment.RData")

mxtvar <- apply(sep_boot_sim$rates, c(1, 2), var)

mxtvar <-mxtvar[reference_age-min(ages_fit),]

dt_all_rates <- data.frame(
  true_rates=true_rates,
  conf1_ifreal= our_rates+simulated_variance_if_real,
  conf2_ifreal= our_rates-simulated_variance_if_real,
  credibility_rates = c(zeros,predicted_rates),
  separate_rates = c(zeros,muxt_hat_predicted_separate),
  conf1_cred = c(zeros,predicted_rates+sqrt(msep)),
  conf2_cred = c(zeros,predicted_rates-sqrt(msep)),
  conf1_sep = c(zeros,muxt_hat_predicted_separate+sqrt(mxtvar)),
  conf2_sep = c(zeros,muxt_hat_predicted_separate-sqrt(mxtvar))
)

# dt_all_rates <- apply(dt_all_rates,2,log)

# library(tidyr)
# library(dplyr)
temporary_data_2 <- dt_all_rates %>% select(true_rates) %>%
  mutate(calendar_time=0:(length(true_rates)-1))



# Assume your dataframe is called dt_all_rates
dt_long <- dt_all_rates %>%  as.data.frame()%>%
  mutate(calendar_time=0:(length(true_rates)-1)) %>%
  select(-true_rates) %>%
  as.data.frame()%>%
  pivot_longer(cols = c(credibility_rates, separate_rates,
                        conf1_ifreal, conf2_ifreal,
                        conf1_cred, conf2_cred, conf1_sep, conf2_sep #,true_rates

                        ),
               names_to = "quantity",
               values_to = "values")

dt_long[['ltype']] <- "solid"
dt_long$ltype[grepl("conf",dt_long$quantity)] <- "dotted"


ggplot(dt_long, aes(x = calendar_time, y = values, color = quantity)) +
  geom_line(aes(linetype = ltype),
            size = 1) +
  geom_point(data=temporary_data_2 %>% rename(values=true_rates) %>% mutate(quantity="true_rates"),
             aes(x=calendar_time,
                 y=values),
             color="gray")+
  scale_linetype_manual(values = c("solid" = "solid", "dotted" = "dotted")) +
  scale_color_manual(values = c(#"true_rates" = "black",
                                "conf1_ifreal" = "black",
                                "conf2_ifreal" = "black",
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

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","msep_65_g2.pdf"),
       width = 8,
       height= 5)

## Reference group 3 ----
reference_age=65
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

insample_rates <- zed*fitted_muxt+(1-zed)*thetaxi*fitted_muxt

our_rates <- c(insample_rates,
               predicted_rates)

tmp_exp_future <- data_pp$list_of_extra_exposures[[reference_group-1]]$Ext[reference_age-min(ages_fit),(ll-forecasting_horizon):ll]

{
  set.seed(1)
  simulation_of_observed_bundles <- sapply(c(tmp_exp,tmp_exp_future)*our_rates, function(x) rpois(1000, x))

  tmp_exp_matrix <- matrix(rep(c(tmp_exp,tmp_exp_future),
                               1000),
                           nrow = 1000,
                           byrow = T)


  simulation_of_observed_bundles <- simulation_of_observed_bundles/tmp_exp_matrix

  simulated_variance_if_real <- apply(simulation_of_observed_bundles, MARGIN = 2, sd)


}





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


# separate_bootstrap_model <- bootstrap(sep_model_fit, nBoot = 500, type = "residual")

# sep_boot_sim<- simulate(separate_bootstrap_model, h = 5)

# save(sep_boot_sim, file = "~/credibility-for-mortality-2024/output/bootg3.RData")

load("~/credibility-for-mortality-2024/output/bootg3.RData")

# save.image("~/credibility-for-mortality-2024/output/my_environmentg3.RData")
#
# load("~/credibility-for-mortality-2024/output/my_environmentg3.RData")

mxtvar <- apply(sep_boot_sim$rates, c(1, 2), var)

mxtvar <-mxtvar[reference_age-min(ages_fit),]

dt_all_rates <- data.frame(
  true_rates=true_rates,
  conf1_ifreal= our_rates+simulated_variance_if_real,
  conf2_ifreal= our_rates-simulated_variance_if_real,
  credibility_rates = c(zeros,predicted_rates),
  separate_rates = c(zeros,muxt_hat_predicted_separate),
  conf1_cred = c(zeros,predicted_rates+sqrt(msep)),
  conf2_cred = c(zeros,predicted_rates-sqrt(msep)),
  conf1_sep = c(zeros,muxt_hat_predicted_separate+sqrt(mxtvar)),
  conf2_sep = c(zeros,muxt_hat_predicted_separate-sqrt(mxtvar))
)

# dt_all_rates <- apply(dt_all_rates,2,log)




# library(tidyr)
# library(dplyr)

temporary_data_2 <- dt_all_rates %>% select(true_rates) %>%
  mutate(calendar_time=0:(length(true_rates)-1))



# Assume your dataframe is called dt_all_rates
dt_long <- dt_all_rates %>%  as.data.frame()%>%
  mutate(calendar_time=0:(length(true_rates)-1)) %>%
  select(-true_rates) %>%
  as.data.frame()%>%
  pivot_longer(cols = c(credibility_rates, separate_rates,
                        conf1_ifreal, conf2_ifreal,
                        conf1_cred, conf2_cred, conf1_sep, conf2_sep #,true_rates

  ),
  names_to = "quantity",
  values_to = "values")

dt_long[['ltype']] <- "solid"
dt_long$ltype[grepl("conf",dt_long$quantity)] <- "dotted"


ggplot(dt_long, aes(x = calendar_time, y = values, color = quantity)) +
  geom_line(aes(linetype = ltype),
            size = 1) +
  geom_point(data=temporary_data_2 %>% rename(values=true_rates) %>% mutate(quantity="true_rates"),
             aes(x=calendar_time,
                 y=values),
             color="gray")+
  scale_linetype_manual(values = c("solid" = "solid", "dotted" = "dotted")) +
  scale_color_manual(values = c(#"true_rates" = "black",
    "conf1_ifreal" = "black",
    "conf2_ifreal" = "black",
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

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","msep_65_g3.pdf"),
       width = 8,
       height= 5)

## Reference group 1 ----
reference_age=65
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

insample_rates <- zed*fitted_muxt+(1-zed)*thetaxi*fitted_muxt

our_rates <- c(insample_rates,
               predicted_rates)

tmp_exp_future <- data_pp$datahat$Ext[reference_age-min(ages_fit),(ll-forecasting_horizon):ll]

{
  set.seed(1)
  simulation_of_observed_bundles <- sapply(c(tmp_exp,tmp_exp_future)*our_rates, function(x) rpois(1000, x))

  tmp_exp_matrix <- matrix(rep(c(tmp_exp,tmp_exp_future),
                               1000),
                           nrow = 1000,
                           byrow = T)


  simulation_of_observed_bundles <- simulation_of_observed_bundles/tmp_exp_matrix

  simulated_variance_if_real <- apply(simulation_of_observed_bundles, MARGIN = 2, sd)


}


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


# separate_bootstrap_model <- bootstrap(sep_model_fit, nBoot = 500, type = "residual")

# sep_boot_sim<- simulate(separate_bootstrap_model, h = 5)

# save(sep_boot_sim,file= "~/credibility-for-mortality-2024/output/bootg1.RData")
load("~/credibility-for-mortality-2024/output/bootg1.RData")

mxtvar <- apply(sep_boot_sim$rates, c(1, 2), var)

mxtvar <-mxtvar[reference_age-min(ages_fit),]

dt_all_rates <- data.frame(
  true_rates=true_rates,
  conf1_ifreal= our_rates+simulated_variance_if_real,
  conf2_ifreal= our_rates-simulated_variance_if_real,
  credibility_rates = c(zeros,predicted_rates),
  separate_rates = c(zeros,muxt_hat_predicted_separate),
  conf1_cred = c(zeros,predicted_rates+sqrt(msep)),
  conf2_cred = c(zeros,predicted_rates-sqrt(msep)),
  conf1_sep = c(zeros,muxt_hat_predicted_separate+sqrt(mxtvar)),
  conf2_sep = c(zeros,muxt_hat_predicted_separate-sqrt(mxtvar))
)

# dt_all_rates <- apply(dt_all_rates,2,log)

dt_all_rates[['calendar_time']] <- 0:(length(true_rates)-1)

# library(tidyr)
# library(dplyr)

temporary_data_2 <- dt_all_rates %>% select(true_rates) %>%
  mutate(calendar_time=0:(length(true_rates)-1))



# Assume your dataframe is called dt_all_rates
dt_long <- dt_all_rates %>%  as.data.frame()%>%
  mutate(calendar_time=0:(length(true_rates)-1)) %>%
  select(-true_rates) %>%
  as.data.frame()%>%
  pivot_longer(cols = c(credibility_rates, separate_rates,
                        conf1_ifreal, conf2_ifreal,
                        conf1_cred, conf2_cred, conf1_sep, conf2_sep #,true_rates

  ),
  names_to = "quantity",
  values_to = "values")

dt_long[['ltype']] <- "solid"
dt_long$ltype[grepl("conf",dt_long$quantity)] <- "dotted"


ggplot(dt_long, aes(x = calendar_time, y = values, color = quantity)) +
  geom_line(aes(linetype = ltype),
            size = 1) +
  geom_point(data=temporary_data_2 %>% rename(values=true_rates) %>% mutate(quantity="true_rates"),
             aes(x=calendar_time,
                 y=values),
             color="gray")+
  scale_linetype_manual(values = c("solid" = "solid", "dotted" = "dotted")) +
  scale_color_manual(values = c(#"true_rates" = "black",
    "conf1_ifreal" = "black",
    "conf2_ifreal" = "black",
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

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","msep_65_g1.pdf"),
       width = 8,
       height= 5)

# Section 4.5 ----
## Assess the model performance in age breaket

ages_fit <- 15:85
years_fit_basic <- 118:140
N_groups=3
forecasting_horizon=1

ageb_width <- 5
out <- NULL

for (seed_ix in 1:3) {

  dt_0 <- data_generator_hmd_lt(
    seed_input = seed_ix,
    hmd_username="",
    hmd_password="",
    exposure_sup = 94500,
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
      answer_sp=NULL
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

        performance_superpop <- model_assessment_superpop(
          model_fit_and_prediction = total_model,
          N_groups = N_groups,
          years_fit = years_fit,
          forecast = forecasting_horizon,
          subset_rates = subset_rates,
          bias=min(ages_fit)
        )

        performance_sp_df <- do.call(rbind,
                                        lapply(performance_superpop, data.frame, stringsAsFactors = FALSE))
        performance_sp_df[['seed']] <- seed_ix
        performance_sp_df[['model']] <- 'superpopulation'
        performance_sp_df[['predictor']] <- rownames(performance_sp_df)
        performance_sp_df[['years_ix']] <- years_ix
        performance_sp_df[['forecasting_horizon']] <- forecasting_horizon
        performance_sp_df[['age_breakets']] <- interval_name
        rownames(performance_sp_df)<-NULL

        answer_sp<-rbind(answer_sp,performance_sp_df)

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
        answer_sp,
        answer_credibility,
        answer_mle,
        answer_separate
      )


    }

}


data.table::fwrite(out,sprintf("~/credibility-for-mortality-2024\\output\\results_age_breakets_smooth_%s.txt",format(Sys.time(), "%Y_%m_%d_%H_%M")))


## Table 1 ----

ages_fit <- 15:85
years_fit_basic <- 118:140
N_groups=3
forecasting_horizon=1

ageb_width <- 5
out <- NULL

for (seed_ix in 1:3) {

  dt_0 <- data_generator_hmd_lt(
    seed_input = seed_ix,
    hmd_username="",
    hmd_password="",
    exposure_sup = 94500,
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

    for(model_ix in c("lc","apc","rh")){

      out <- rbind(out,
                   data.frame(
                     model=model_ix,
                     years_index=years_ix,
                     bic=BIC(mortality_models_fit[[model_ix]])
                   ))

    }
}

}


data.table::fwrite(out,"~/credibility-for-mortality-2024\\output\\out_insample_bic.txt")


# Import ----
library(data.table)
library(xtable)
library(dplyr)

source(
  "~/credibility-for-mortality-2024/data_generators.R"
)

source(
  "~/credibility-for-mortality-2024/data_preprocessing.R"
)

source(
  "~/credibility-for-mortality-2024/utils_credibility.R"
)

# Simulation study -----

## Generate the data, pre-process and fit the global mortality model ----

ages_fit <- 15:85
years_fit<- 118:140
N_groups=3


### data generation
dt_0 <- data_generator_hmd_lt(
  seed_input = 1,
  hmd_username="",
  hmd_password="",
  exposure_sup = 100000,
  exposure_sub = c(5000, 500)
)


### data pre-processing
data_pp <- data_preprocessing_scenario_2(
  data = dt_0,
  ages_fit = ages_fit,
  years_fit=years_fit
  
)

### global mortality model fit
mortality_models_fit <- fit_mortality_models(data_pp, years_fit-min(years_fit), ages_fit,separate_exposures =F)


## Figure 5 ----
## Numerical application in age-groups 
dt_ab <- fread('~/GitHub/credibility-for-mortality-2024/output/results_age_breakets_smooth_2025_06_16_14_56.txt')

dt_out=dt_ab %>%
  filter(model != "total") %>%
  group_by(seed,model, predictor, age_breakets) %>%
  summarise(error=mean(error),
            oos_deviance=mean(oos_deviance)) %>%
  group_by(model, predictor, age_breakets) %>%
  reframe(error_resp=mean(error),
          error_sd=sd(error),
          oos_deviance_resp=mean(oos_deviance),
          oos_deviance_sd=sd(oos_deviance))

# dt_out %>%
#   filter(predictor  =="apc",
#          age_breakets=="16-25")

ggplot(dt_out %>%
         filter(predictor  =="lc"), aes(x = age_breakets, y = error_resp, color = model, group = model)) +
  geom_line(linewidth=2) +    # Add lines
  geom_point() +   # Add points
  geom_errorbar(aes(ymin = error_resp-error_sd ,
                    ymax = error_resp+error_sd ),width = 0.2,alpha=.3)+   # Add points
  theme_minimal(base_size = 28) +
  ylim(0,1)+
  labs(x = "", y = "", color = "")+
  scale_color_discrete(labels = c("A", "B","C","D")) +
  theme(      legend.position = "inside",                          # must now be "inside"
              # legend.position.inside = c(0.85, 0.85),              # numeric position
              legend.justification = c("right", "top"),
              legend.key.height = unit(1.5, "lines"),
              legend.background = element_blank(),   # removes legend box
              legend.key = element_blank(),
              axis.text.x = element_text(angle = 90))

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_mare_age_classes.pdf"),
       width = 8,
       height= 5)

ggplot(dt_out %>%
         filter(predictor  =="lc"), aes(x = age_breakets, y = oos_deviance_resp, color = model, group = model)) +
  geom_line(size=2) +    # Add lines
  geom_point() +   # Add points
  geom_errorbar(aes(ymin = oos_deviance_resp-oos_deviance_sd,
                    ymax = oos_deviance_resp+oos_deviance_sd),width = 0.2,alpha=.3)+
  theme_minimal(base_size = 28) +
  ylim(0,30)+
  labs(x = "", y = "", color = "")+
  scale_color_discrete(labels = c("A", "B","C","D")) +
  theme(      legend.position = "inside",                          # must now be "inside"
              # legend.position.inside = c(0.85, 0.85),              # numeric position
              legend.justification = c("left", "top"),
              legend.key.height = unit(1.5, "lines"),
              legend.background = element_blank(),   # removes legend box
              legend.key = element_blank(),
              axis.text.x = element_text(angle = 90))

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_deviance_age_classes.pdf"),
       width = 8,
       height= 5)


# Same plots, now by group ----
# ~/credibility-for-mortality-2024/output/results_age_breakets_smooth_2025_06_04_10_17.txt
dt_ab_pgroup <- fread('~/GitHub/credibility-for-mortality-2024/output/results_age_breakets_smooth_2025_06_16_14_56.txt')

## Group 0 ----

dt_out=dt_ab_pgroup %>%
  filter(model != "superpopulation") %>%
  group_by(seed,model, predictor, age_breakets) %>%
  summarise(error=mean(error_1),
            oos_deviance=mean(oos_deviance_1)) %>%
  group_by(model, predictor, age_breakets) %>%
  reframe(error_resp=mean(error),
          error_sd=sd(error),
          oos_deviance_resp=mean(oos_deviance),
          oos_deviance_sd=sd(oos_deviance))


ggplot(dt_out %>%
         filter(predictor  =="lc"), aes(x = age_breakets, y = error_resp, color = model, group = model)) +
  geom_line(size=2) +    # Add lines
  geom_point() +   # Add points
  geom_errorbar(aes(ymin = error_resp-error_sd ,
                    ymax = error_resp+error_sd ),width = 0.2,size=2,alpha=.3)+
  theme_minimal(base_size = 28) +
  ylim(0,1)+
  labs(x = "", y = "", color = "")+
  scale_color_discrete(labels = c("A", "B","C","D")) +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "inside",                          # must now be "inside"
    # legend.position.inside = c(0.85, 0.85),              # numeric position
    legend.justification = c("right", "top"),
    legend.key.height = unit(1.5, "lines"),
    legend.background = element_blank(),   # removes legend box
    legend.key = element_blank())

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_mare_age_classes_group_1.pdf"),
       width = 8,
       height= 7)

ggplot(dt_out %>%
         filter(predictor  =="lc"), aes(x = age_breakets, y = oos_deviance_resp, color = model, group = model)) +
  geom_line(size=2) +    # Add lines
  geom_point() +   # Add points
  geom_errorbar(aes(ymin = oos_deviance_resp-oos_deviance_sd ,
                    ymax = oos_deviance_resp+oos_deviance_sd ),width = 0.2,size=2,alpha=.3)+
  theme_minimal(base_size = 28) +
  ylim(0,30)+
  labs(x = "", y = "", color = "")+
  scale_color_discrete(labels = c("A", "B","C","D")) +
  theme(      legend.position = "inside",                          # must now be "inside"
              # legend.position.inside = c(0.85, 0.85),              # numeric position
              legend.justification = c("left", "top"),
              legend.key.height = unit(1.5, "lines"),
              legend.background = element_blank(),   # removes legend box
              legend.key = element_blank(),
              axis.text.x = element_text(angle = 90))

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_deviance_age_classes_group_1.pdf"),
       width = 8,
       height= 7)

## Group 1 ----

dt_out=dt_ab_pgroup %>%
  filter(model != "superpopulation") %>%
  group_by(seed,model, predictor, age_breakets) %>%
  summarise(error=mean(error_2,na.rm=T),
            oos_deviance=mean(oos_deviance_2,na.rm=T)) %>%
  group_by(model, predictor, age_breakets) %>%
  reframe(error_resp=mean(error),
          error_sd=sd(error),
          oos_deviance_resp=mean(oos_deviance),
          oos_deviance_sd=sd(oos_deviance))


ggplot(dt_out %>%
         filter(predictor  =="lc"), aes(x = age_breakets, y = error_resp, color = model, group = model)) +
  geom_line(size=2) +    # Add lines
  geom_point() +   # Add points
  geom_errorbar(aes(ymin = error_resp-error_sd ,
                    ymax = error_resp+error_sd ),width = 0.2,size=2,alpha=.3)+
  theme_minimal(base_size = 28) +
  ylim(0,1)+
  labs(x = "", y = "", color = "")+
  scale_color_discrete(labels = c("A", "B","C","D")) +
  theme(      legend.position = "inside",                          # must now be "inside"
              # legend.position.inside = c(0.85, 0.85),              # numeric position
              legend.justification = c("right", "top"),
              legend.key.height = unit(1.5, "lines"),
              legend.background = element_blank(),   # removes legend box
              legend.key = element_blank(),
              axis.text.x = element_text(angle = 90))

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_mare_age_classes_group_2.pdf"),
       width = 8,
       height= 7)


ggplot(dt_out %>%
         filter(predictor  =="lc"), aes(x = age_breakets, y = oos_deviance_resp, color = model, group = model)) +
  geom_line(size=2) +    # Add lines
  geom_point() +   # Add points
  geom_errorbar(aes(ymin = oos_deviance_resp-oos_deviance_sd ,
                    ymax = oos_deviance_resp+oos_deviance_sd ),width = 0.2,size=2,alpha=.3)+
  theme_minimal(base_size = 28) +
  ylim(0,30)+
  labs(x = "", y = "", color = "")+
  scale_color_discrete(labels = c("A", "B","C","D")) +
  theme(      legend.position = "inside",                          # must now be "inside"
              # legend.position.inside = c(0.85, 0.85),              # numeric position
              legend.justification = c("left", "top"),
              legend.key.height = unit(1.5, "lines"),
              legend.background = element_blank(),   # removes legend box
              legend.key = element_blank(),
              axis.text.x = element_text(angle = 90))

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_deviance_age_classes_group_2.pdf"),
       width = 8,
       height= 7)

## Group 2 ----


dt_out=dt_ab_pgroup %>%
  filter(model != "superpopulation") %>%
  group_by(seed,model, predictor, age_breakets) %>%
  summarise(error=mean(error_3,na.rm=T),
            oos_deviance=mean(oos_deviance_3,na.rm=T)) %>%
  group_by(model, predictor, age_breakets) %>%
  reframe(error_resp=mean(error,na.rm=T),
          error_sd=sd(error,na.rm=T),
          oos_deviance_resp=mean(oos_deviance,na.rm=T),
          oos_deviance_sd=sd(oos_deviance,na.rm=T))


ggplot(dt_out %>%
         filter(predictor  =="lc"), aes(x = age_breakets, y = error_resp, color = model, group = model)) +
  geom_line(size=2) +    # Add lines
  geom_point() +   # Add points
  geom_errorbar(aes(ymin = error_resp-error_sd ,
                    ymax = error_resp+error_sd ),width = 0.2,size=2,alpha=.3)+
  theme_minimal(base_size = 28) +
  labs(x = "", y = "", color = "")+
  ylim(0,1)+
  scale_color_discrete(labels = c("A", "B","C","D")) +
  theme(      legend.position = "inside",                          # must now be "inside"
              # legend.position.inside = c(0.85, 0.85),              # numeric position
              legend.justification = c("right", "top"),
              legend.key.height = unit(1.5, "lines"),
              legend.background = element_blank(),   # removes legend box
              legend.key = element_blank(),
              axis.text.x = element_text(angle = 90))

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_mare_age_classes_group_3.pdf"),
       width = 8,
       height= 7)


ggplot(dt_out %>%
         filter(predictor  =="lc") , aes(x = age_breakets, y = oos_deviance_resp, color = model, group = model)) +
  geom_line(size=2) +    # Add lines
  geom_point() +   # Add points
  geom_errorbar(aes(ymin = oos_deviance_resp-oos_deviance_sd,
                    ymax = oos_deviance_resp+oos_deviance_sd),width = 0.2,size=2,alpha=.3)+
  theme_minimal(base_size = 28) +
  ylim(0,30)+
  labs(x = "", y = "", color = "")+
  scale_color_discrete(labels = c("A", "B","C","D")) +
  theme(      legend.position = "inside",                          # must now be "inside"
              # legend.position.inside = c(0.85, 0.85),              # numeric position
              legend.justification = c("left", "top"),
              legend.key.height = unit(1.5, "lines"),
              legend.background = element_blank(),   # removes legend box
              legend.key = element_blank(),
              axis.text.x = element_text(angle = 90))

ggsave(paste0("~/credibility-for-mortality-2024\\output\\","lc","_deviance_age_classes_group_3.pdf"),
       width = 8,
       height= 7)

# remark 4 ----
# Separate models convergence

library(data.table)


folder_path <- "~/credibility-for-mortality-2024/output/"

file_list <- list.files(path = folder_path, pattern = "^failed_fit_.*\\.txt$", full.names = TRUE)

# Read all files into a single data.table
combined_data <- rbindlist(lapply(file_list, fread), fill = TRUE, idcol = "Source_File")

summarized_data <- combined_data[, .(Total_Failed = sum(failed_fit, na.rm = TRUE)), by = .(group,model)]

round(summarized_data[['Total_Failed']]/(3*6),2)









