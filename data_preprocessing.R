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
  
  if (scenario == 2) {
    out <- data_preprocessing_scenario_2(
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



data_preprocessing_scenario_2 <- function(data,
                                          N_groups = 3,
                                          bias=0,
                                          years_fit=110:150,
                                          ages_fit = 50:65) {
  
# browser()
 tmp_occ <-  data %>%
    tidyr::pivot_wider(id_cols=Age,names_from = Period, values_from = Dx) %>%
   as.matrix()
 rownames(tmp_occ) <- tmp_occ[,"Age"]
 tmp_occ <- tmp_occ[,-1]
 tmp_occ <- tmp_occ[as.character(ages_fit),years_fit+bias]
 
 
 tmp_exp <-  data %>%
   tidyr::pivot_wider(id_cols=Age,names_from = Period, values_from = Ext) %>%
   as.matrix()
 rownames(tmp_exp) <- tmp_exp[,"Age"]
 tmp_exp <- tmp_exp[,-1]
 tmp_exp <- tmp_exp[as.character(ages_fit),years_fit+bias]
 colnames(tmp_occ) <-colnames(tmp_exp) <- years_fit-min(years_fit)
 
 tmp_occ1 <-  data %>%
   tidyr::pivot_wider(id_cols=Age,names_from = Period, values_from = Dx1) %>%
   as.matrix()
 rownames(tmp_occ1) <- tmp_occ1[,"Age"]
 tmp_occ1 <- tmp_occ1[,-1]
 tmp_occ1 <- tmp_occ1[as.character(ages_fit),years_fit+bias]
 
 
 tmp_exp1 <-  data %>%
   tidyr::pivot_wider(id_cols=Age,names_from = Period, values_from = Ext1) %>%
   as.matrix()
 rownames(tmp_exp1) <- tmp_exp1[,"Age"]
 tmp_exp1 <- tmp_exp1[,-1]
 tmp_exp1 <- tmp_exp1[as.character(ages_fit),years_fit+bias]
 colnames(tmp_occ1) <-colnames(tmp_exp1) <- years_fit-min(years_fit)
 
 tmp_occ2 <-  data %>%
   tidyr::pivot_wider(id_cols=Age,names_from = Period, values_from = Dx2) %>%
   as.matrix()
 rownames(tmp_occ2) <- tmp_occ2[,"Age"]
 tmp_occ2 <- tmp_occ2[,-1]
 tmp_occ2 <- tmp_occ2[as.character(ages_fit),years_fit+bias]
 
 
 tmp_exp2 <-  data %>%
   tidyr::pivot_wider(id_cols=Age,names_from = Period, values_from = Ext2) %>%
   as.matrix()
 rownames(tmp_exp2) <- tmp_exp2[,"Age"]
 tmp_exp2 <- tmp_exp2[,-1]
 tmp_exp2 <- tmp_exp2[as.character(ages_fit),years_fit+bias]
 colnames(tmp_occ2) <-colnames(tmp_exp2) <- years_fit-min(years_fit)
  
 datahat <- structure(
   list(
     Dxt = tmp_occ,#+tmp_occ1+tmp_occ2,
     Ext = tmp_exp,#+tmp_exp1+tmp_exp2,
     ages = ages_fit,
     years = years_fit-min(years_fit),
     type = 'central',
     series = 'female',
     label = 'total'
   ),
   class = "StMoMoData"
 )
 
 datahat_tot <- structure(
   list(
     Dxt = tmp_occ+tmp_occ1+tmp_occ2,
     Ext = tmp_exp+tmp_exp1+tmp_exp2,
     ages = ages_fit,
     years = years_fit-min(years_fit),
     type = 'central',
     series = 'female',
     label = 'total'
   ),
   class = "StMoMoData"
 )
  
  l <- list()
  
  list_of_extra_exposures <- list()
 
  
  list_of_extra_exposures[[1]] <- list(Dxt=tmp_occ1,
                                       Ext= tmp_exp1)
  

  
  list_of_extra_exposures[[2]] <- list(Dxt=tmp_occ2,
                                       Ext= tmp_exp2)
  
  
  l[['datahat']] <- datahat
  l[['datahat_tot']] <- datahat_tot
  l[['list_of_extra_exposures']] <- list_of_extra_exposures
  
  return(l)
  
}

