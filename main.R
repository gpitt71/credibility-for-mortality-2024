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

ages_fit <- 50:80
N_groups=3

dt_0 <- data_generator_1(seed = 1)

data_pp <- data_preprocessing_scenario_1(
  data = dt_0,
  N_groups = N_groups,
  ages_fit = ages_fit,
  years_fit = 0:29
)

mortality_models_fit <- fit_mortality_models(data_pp, 0:29, ages_fit)

credibility_model <- fit_and_predict_credibility_models(
  data = dt_0,
  data_pp =  data_pp,
  scenario = 1,
  N_groups = N_groups,
  mortality_models_fit = mortality_models_fit,
  years_fit = 0:29,
  ages_fit = ages_fit,
  forecasting_horizon = 1
)


## LC: thetas 75----

thetas_plotter(
  credibility_model,
  subgroup=1,
  ages_fit,
  years_fit,
  mortality_models_fit,
  predictor = "lc",
  chosen_age = 75
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
  chosen_age = 75
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
  chosen_age = 75
)

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_simulation_cs_3.pdf"),
       width = 8,
       height= 5)

## LC: thetas ----

all_thetas_plotter(credibility_model, N_groups, ages_fit, predictor = "lc") +
  scale_color_manual(
    values = c("C1" = "#a71429","C2" = "#4169E1", "C3" = "#006400"),
    labels = c(expression(hat(Theta)[x]^1), expression(hat(Theta)[x]^2), expression(hat(Theta)[x]^3))
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
    labels = c(expression(1 - Z[x]^1), expression(1 - Z[x]^2), expression(1 - Z[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","lc","_simulation_zs.pdf"),
       width = 8,
       height= 5)

## APC: weights ----

p <- weights_plotter (credibility_model, N_groups, ages_fit, predictor = "apc")

p +
  scale_color_manual(
    values = c("Z_1" = "#a71429","Z_2" = "#4169E1", "Z_3" = "#006400"),
    labels = c(expression(1 - Z[x]^1), expression(1 - Z[x]^2), expression(1 - Z[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","apc","_simulation_zs.pdf"),
       width = 8,
       height= 5)

## RH: weights ----

p <- weights_plotter (credibility_model, N_groups, ages_fit, predictor = "rh")

p +
  scale_color_manual(
    values = c("Z_1" = "#a71429","Z_2" = "#4169E1", "Z_3" = "#006400"),
    labels = c(expression(1 - Z[x]^1), expression(1 - Z[x]^2), expression(1 - Z[x]^3))
  )

ggsave(paste0("C:\\Users\\pwt887\\Documents\\GitHub\\credibility-for-mortality-2024\\output\\","rh","_simulation_zs.pdf"),
       width = 8,
       height= 5)





# Some parameters ----
ages_fit <- 50:80
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







