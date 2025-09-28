################################################################################
# IMPORTANT: THE MORTALITY DATA provided are simulated FROM MODEL 4 IN THE
# PAPER (SB-DLNM WITH TIME-SERIES DESIGN) due to confidentiality restrictions of 
# the original data set. Original data may be requested to the data owner (The 
# Barcelona Public Health Agency) who should share them in similar terms than 
# those applying to this study. The supplied mortality, with exactly the same 
# structure than the original data set, allows reproducing the code provided so 
# that interested readers may run it as an example of use. However, WinBUGS 
# results from our original analyses are additionally supplied (input/result_paper
# folder) so that readers should be able to reproduce exactly our analyses after 
# those WinBUGS calls (model summaries and plots).
################################################################################

################################################################################
# In this R project we implement Bayesian and Spatial Bayesian Distributed 
# Lag Non-Linear Models (B-DLNM and SB-DLNM) for the case Study of short-term 
# associations between temperature and mortality in the city of Barcelona.
################################################################################

################################################################################
# CODE 1: DATA PREPARATION
# Data preparation: Transforming time-series temperature and mortality data sets
# for B-DLNMs and SB-DLNMs.
################################################################################

# Load libraries

library(sf)
library(spdep)
library(lubridate)
library(dlnm)
library(gnm)
library(splines)
library(rstan)
library(cmdstanr)
library(tidyverse)
library(foreign)
library(shinystan)
# Load data

load("input/daily_data.RData")
head(data)

# Set variables defining the dlnm model
load("output/dlnm_configuration.RData")
dlnm_var  

# Set variables for trend and seasonality

# Subset data to only summer months of 2007 to 2016
data <- subset(data, month(date) %in% 6:9)
data <- subset(data, year(date) %in% 2007:2016)

data$year <- factor(lubridate::year(data$date))
head(data)
# Create crossbasis for each region

# Ensure that the data is ordered by region to maintain alignment between 
# crossbasis and regions

if(is.unsorted(data$region)) {
  stop("data needs to be ordered by region for the next loop")}

list_X <- vector("list", dlnm_var$n_reg)

for(i_reg in 1:dlnm_var$n_reg) {
  
      ## this is all the same as in 01_data_preparation.R
      
      ## >> needed to make this sprintf because `region` is a string
      temp <- subset(data, region == sprintf("%02i", i_reg), 
                     select = c("temp", paste0("lag", 1:dlnm_var$max_lag)))
      
      temp_knots <- quantile(temp$temp, dlnm_var$var_prc, na.rm = TRUE)
      temp_boundary <- range(temp, na.rm = TRUE)
      
      cb <- crossbasis(temp,
                       argvar = list(fun = dlnm_var$var_fun,
                                     knots = temp_knots,
                                     Boundary.knots = temp_boundary),
                       arglag = list(fun = "ns",
                                     knots = logknots(dlnm_var$max_lag, 
                                                      dlnm_var$lagnk),
                                     intercept = TRUE))

    
  ##
  list_X[[i_reg]] <- cb
}

# PREPARE DATA FOR THE CASE-CROSSOVER DESIGN

# Create strata for the case-crossover
# (neighborhood - year - month - day of week)
data$strata <- paste(data$region, 
                     year(data$date), 
                     formatC(month(data$date), width = 2, flag = "0"),
                     wday(data$date, week_start = 1),
                     sep = ":")

# make them both lists
data_l <- split(data, f = data$region)

#' ===========================================================================
#' Normal Cond and Poisson

library(gnm)

# just try 1 area that seems like it has a strong relationship
# region 39
df_local <- data_l[[39]]
cb_local <- list_X[[39]]

# keep just the non-empty strata
df_local_agg <- df_local %>%
  group_by(strata) %>%
  summarize(
    .groups = 'keep',
    sum_mort = sum(mort)
  ) %>% 
  mutate(keep = ifelse(sum_mort > 0, 1, 0))

df_local <- left_join(df_local, df_local_agg)
table(df_local$keep, useNA = 'always')

# quasi-poisson
m1.std <- gnm(mort ~ cb_local + factor(strata), 
              data = df_local,
              family = quasipoisson,
              subset = keep == 1)

cp1 <- crosspred(cb_local, m1.std)
plot(cp1, 'overall')

# but are the same as the base model
load("input/result_paper/final_simsmatrix_model1_independent_casecrossover.RData")
for(i_reg in 1:73) {
  beta_reg1 <- winbugs_res[,grepl(paste0("^beta\\[", i_reg,","), 
                                  colnames(winbugs_res))]
  
  xx <- cbind("paper" = apply(beta_reg1, 2, median), 
        "freq-P" = coef(m1.std)[2:13] )
  
  if(sum(apply(xx, 1, function(a) abs(a[1] - a[2])^2)) < 1 ) stop()
  
  print(xx)
}
rm(winbugs_res); gc();

load("input/result_paper/final_simsmatrix_model2_independent_timseseries.RData")
for(i_reg in 1:73) {
  beta_reg1 <- winbugs_res[,grepl(paste0("^beta\\[", i_reg,","), 
                                  colnames(winbugs_res))]
  xx <- cbind("paper" = apply(beta_reg1, 2, median), 
              "freq-P" = coef(m1.std)[2:13] )
  if(sum(apply(xx, 1, function(a) abs(a[1] - a[2])^2)) < 1 ) stop()
  print(xx)
}
rm(winbugs_res); gc();


