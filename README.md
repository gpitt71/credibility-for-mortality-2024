# Forecasting sub-population mortality using credibility theory

This repository contains the code to the results from the manuscript *Forecasting sub-population mortality using credibility theory*.

All the manuscript results can be replicated running the `main.R` script.

In order to run the `main.R` script you will need to create an account in the [Human Mortality Database](https://www.mortality.org/) website.

Inside the `main.R` you might need to:

* Pass the username and password that you created as arguments of class `character` to the `data_generator_hmd_lt` function.
* Change the directories pattern to a pattern that works on your computer.


The content of the repository is:

- `main.R` can be used to reproduce our models results. The sections of the code refer to the sections of the manuscript where the code has been used..
- `data_generators.R` contains the function to generate the data.
- `data_preprocessing.R` contains the function to pre-process the data. 
- `utils_credibility.R` contains the functions that we use behind the curtains to perform the computations.
-  The `output` folder where the results will be stored. In the output folder we included three pre-saved R environment that will be loaded while executing the script to compare the mean squared error of prediction of our model against model C. The three environments include the results of pre-simulated prediction error for model C in the different groups using the procedure outlined in the main paper.

# R system and packages used

```
R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8    LC_MONETARY=Danish_Denmark.utf8
[4] LC_NUMERIC=C                    LC_TIME=Danish_Denmark.utf8    

time zone: Europe/Copenhagen
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] rpart_4.1.24      data.table_1.16.4 tidyr_1.3.1       ggplot2_3.5.2     StMoMo_0.4.1.9000 forecast_8.23.0  
[7] gnm_1.1-5         dplyr_1.1.4      

loaded via a namespace (and not attached):
 [1] generics_0.1.3     lattice_0.22-6     digest_0.6.37      magrittr_2.0.3     evaluate_1.0.3     grid_4.4.2        
 [7] RColorBrewer_1.1-3 maps_3.4.2.1       fastmap_1.2.0      Matrix_1.7-1       nnet_7.3-19        purrr_1.0.4       
[13] spam_2.11-0        viridisLite_0.4.2  scales_1.4.0       fanplot_4.0.0      cli_3.6.3          rlang_1.1.4       
[19] withr_3.0.2        rootSolve_1.8.2.4  tools_4.4.2        parallel_4.4.2     colorspace_2.1-1   curl_6.0.1        
[25] vctrs_0.6.5        R6_2.6.1           zoo_1.8-12         lifecycle_1.0.4    tseries_0.10-58    relimp_1.0-5      
[31] MASS_7.3-61        pkgconfig_2.0.3    urca_1.3-4         pillar_1.10.2      gtable_0.3.6       glue_1.8.0        
[37] quantmod_0.4.26    Rcpp_1.0.13-1      fields_16.3        xfun_0.52          tibble_3.2.1       lmtest_0.9-40     
[43] tidyselect_1.2.1   rstudioapi_0.17.1  knitr_1.50         farver_2.1.2       qvcalc_1.0.3       htmltools_0.5.8.1 
[49] nlme_3.1-166       rmarkdown_2.29     xts_0.14.1         dotCall64_1.2      timeDate_4041.110  fracdiff_1.5-3    
[55] compiler_4.4.2     quadprog_1.5-8     TTR_0.24.4 
```
