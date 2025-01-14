data_preprocessing <- function(data,
                               scenario = 0,
                               N_groups = 3,
                               ages_fit = 50:80,
                               years_fit = 0:29) {
  if (scenario == 0) {
    out <- data_preprocessing_scenario_0(
      data,
      N_groups = N_groups,
      ages_fit = ages_fit,
      years_fit = years_fit
    )
    
    
  }
  
  if (scenario == 1) {
    out <- data_preprocessing_scenario_1(
      data,
      N_groups = N_groups,
      ages_fit = ages_fit,
      years_fit = years_fit
    )
    
    
  }
  
  
  return(out)
  
}




data_preprocessing_scenario_0 <- function(data,
                                          N_groups = 3,
                                          ages_fit = 50:80,
                                          years_fit = 0:29) {
  l <- list()
  
  list_of_extra_exposures <- list()
  
  superpop_data <- data$sim_out_superp
  subpop_data <- data$sim_out_subp
  
  # super-pop class is set to one
  
  occurrence_tab <- death_table(
    superpop_data$population[superpop_data$population$risk_cls == 1, ],
    ages = c(ages_fit, last(ages_fit) + 1),
    period = c(years_fit, last(years_fit) + 1)
  )
  
  exposure_tab <- exposure_table(
    superpop_data$population[superpop_data$population$risk_cls == 1, ],
    ages = c(ages_fit, last(ages_fit) + 1),
    period = c(years_fit, last(years_fit) + 1)
  )
  
  
  datahat <- structure(
    list(
      Dxt = occurrence_tab,
      Ext = exposure_tab,
      ages = ages_fit,
      years = years_fit,
      type = 'central',
      series = 'female',
      label = 'total'
    ),
    class = "StMoMoData"
  )
  
  
  # Extra exposure data
  
  for (i in 2:N_groups) {
    occurrence_tab <- death_table(
      subpop_data$population[subpop_data$population$risk_cls == i, ],
      ages = c(ages_fit, last(ages_fit) + 1),
      period = c(years_fit, last(years_fit) + 1)
    )
    
    exposure_tab <- exposure_table(
      subpop_data$population[subpop_data$population$risk_cls == i, ],
      ages = c(ages_fit, last(ages_fit) + 1),
      period = c(years_fit, last(years_fit) + 1)
    )
    
    
    assign(paste0("E", i), exposure_tab)
    assign(paste0("D", i), occurrence_tab)
    
    tmp_list <- list()
    
    tmp_list[['Dxt']] <- get(paste0("D", i))
    tmp_list[['Ext']] <- get(paste0("E", i))
    list_of_extra_exposures[[i - 1]] <- tmp_list
    
  }
  
  
  l[['datahat']] <- datahat
  l[['list_of_extra_exposures']] <- list_of_extra_exposures
  
  return(l)
  
}


data_preprocessing_scenario_1 <- function(data,
                                          N_groups = 3,
                                          ages_fit = 50:80,
                                          years_fit = 0:29) {
  l <- list()
  
  list_of_extra_exposures <- list()
  
  pop_data <- data
  
  # super-pop class is set to one
  
  occurrence_tab <- death_table(
    pop_data$population[pop_data$population$risk_cls == 1, ],
    ages = c(ages_fit, last(ages_fit) + 1),
    period = c(years_fit, last(years_fit) + 1)
  )
  
  exposure_tab <- exposure_table(
    pop_data$population[pop_data$population$risk_cls == 1, ],
    ages = c(ages_fit, last(ages_fit) + 1),
    period = c(years_fit, last(years_fit) + 1)
  )
  
  
  datahat <- structure(
    list(
      Dxt = occurrence_tab,
      Ext = exposure_tab,
      ages = ages_fit,
      years = years_fit,
      type = 'central',
      series = 'female',
      label = 'total'
    ),
    class = "StMoMoData"
  )
  
  
  # Extra exposure data
  
  for (i in 2:N_groups) {
    occurrence_tab <- death_table(
      pop_data$population[pop_data$population$risk_cls == i, ],
      ages = c(ages_fit, last(ages_fit) + 1),
      period = c(years_fit, last(years_fit) + 1)
    )
    
    exposure_tab <- exposure_table(
      pop_data$population[pop_data$population$risk_cls == i, ],
      ages = c(ages_fit, last(ages_fit) + 1),
      period = c(years_fit, last(years_fit) + 1)
    )
    
    
    assign(paste0("E", i), exposure_tab)
    assign(paste0("D", i), occurrence_tab)
    
    tmp_list <- list()
    
    tmp_list[['Dxt']] <- get(paste0("D", i))
    tmp_list[['Ext']] <- get(paste0("E", i))
    list_of_extra_exposures[[i - 1]] <- tmp_list
    
  }
  
  
  l[['datahat']] <- datahat
  l[['list_of_extra_exposures']] <- list_of_extra_exposures
  
  return(l)
  
}
