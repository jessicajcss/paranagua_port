
###############################################################################
# Title: Case-Crossover of Port (Un)Berthing and Moored Activities vs Air Quality
# Date: 2025-08-17 [original: C:/Users/air_p/OneDrive - ufpr.br/week_temps/rmc_relatorio/outros projetos/porto/Porto - Rafaela/Porto - graph/scripts/FINAL/aplicado/20250904/03-case-crossoverCoPilot_berthing_moored_baselinedataset_diagnostic_plot.R]
# Last update: 2025-09-05 
# Overview (methods in brief)
# ---------------------------
# We use time-stratified case-crossover designs (Maclure 1991; Janes et al. 2005):
#   • Within each stratum (site × hour-of-day × day-of-week × time-block), we compare
#     minutes/intervals that meet a "case" definition to nearby "control" intervals.
#   • Strata serve as fixed effects, removing all time-invariant confounding within
#     the block (seasonality, slow trends, sensor baselines, etc.).
#   • We form cases using sensor-specific quantile thresholds (e.g., 75th percentile),
#     which respects vertical gradient/scale differences across ALPHA/BETA/GAMMA.
#
# Two target analyses
# -------------------
# A) (UN)BERTHING EPISODES:
#    Case = pollutant > sensor-specific threshold AND within ± 20 min of a berthing event [berthing > 0]
#    Strata = site × hour × DOW × month
#
# B) MOORED OPERATIONS (loading/unloading/handling):
#    Case = pollutant > sensor-specific threshold  AND matching cargo
#    Strata = site × hour × DOW × month
#
# Design choices & diagnostics we implement
# -----------------------------------------
# • Time-stratified matching (Janes/Sheppard/Lumley 2005) to avoid overlap bias.
# • Cap controls per case (e.g., ≤ 4:1) to prevent quasi-separation and variance inflation
#   with extremely imbalanced strata (rule-of-thumb; see also EPV literature).
# • Minimum number of cases (e.g., ≥ 40–50) and of usable strata (≥20) to ensure stable
#   conditional logistic estimates (Harrell 2015; Vittinghoff & McCulloch 2007).
# • Report AIC/BIC for transparency of in-sample fit (keep identical covariates across
#   specs when using AIC/BIC for model ranking).
#
# Key references (general, accessible)
# ------------------------------------
# Maclure M. (1991). The case-crossover design: a method for studying transient effects
#   on the risk of acute events. Am J Epidemiol, 133(2):144–153.
# Navidi W. (1998). Bidirectional case-crossover designs for exposures with time trends.
#   Biometrics, 54(2):596–605.
# Janes H, Sheppard L, Lumley T. (2005). Case-crossover analyses of air pollution exposure
#   data: referent selection strategies and their implications for bias. Epidemiology, 16(6):717–726.
# Mittleman MA. (2005). Optimal referent selection strategies in case-crossover studies.
#   Am J Epidemiol, 161(2):161–168.
# Bateson TF, Schwartz J. (1999). Control for seasonal variation and weather in time-series
#   models of air pollution and mortality. Environ Res, 81(3): 297–309.  (context)
# Vittinghoff E, McCulloch CE. (2007). Relaxing the rule of ten events per variable in
#   logistic and Cox regression. Am J Epidemiol, 165(6):710–718.
# Harrell FE. (2015). Regression Modeling Strategies. Springer.  (EPV & model stability)
###############################################################################
# ===============================================================
# This script runs primary and sensitivity analyses for both models:
#   MODEL A: (Un)Berthing episodes
#   MODEL B: Moored operations
# and outputs a combined results table with specification labels.
# ===============================================================


# Author: Santos-Silva, J. C.
# Last Updated: 05/09/2025
################################################################################



library(lubridate)
library(bazar)
library(rstatix)
library(dplyr)
library(tidyverse)
library(openair)
library(data.table)


################################################################################
################################################################################
########## Previous run: Data INPUT ---
################################################################################
################################################################################


#==== METEO ----
meteo <- fread('./data/geral/meteo.csv') %>%# local time
  mutate(date = dmy_hm(date))

meteo_simepar <- data.table::fread("./data/dados_simepar.csv") %>%
  mutate(date = as_datetime(datahora)) %>%
  select(-datahora)


#==== TRAJ MIXDEPTH ----
trajetorias_10 <- readRDS("output/paranagua_96h_10m_60min.rds")


traj_10 <- trajetorias_10 %>%
  selectByDate(start = "07/09/2024", 
               end ="03/12/2024") %>%
  # Convert to GMT-3 (e.g., America/Sao_Paulo)
  mutate(date_utc = date,
         date = with_tz(date_utc, tzone = "America/Sao_Paulo")) %>%
  subset(hour.inc == 0) %>% 
  select(date, mixdepth)


#==== POLLUTANTS ----
pollutants <- fread('./data/geral/pollutants.csv') %>% # time TZ; pm in ug/m³ & gases in ppb
  mutate(date = dmy_hm(date),
         site = case_when(sensor == "alfa" ~ "ALPHA", 
                          sensor == "gamma" ~ "BETA", # GAMMA and BETA were interchanged according to site position indicated in the paper
                          sensor == "beta" ~ "GAMMA")) %>%
  subset(date >= ymd_hms("2024-09-07 00:00:00") & 
           date <= ymd_hms("2024-12-04 23:59:00")) %>%
  # correction to ug/m³
  mutate(no2 = no2*1.88,
         so2 = so2*2.62)



#==== ATRACACOES ---- 

# DOCKING --
atracacoes1 <- fread('./data/geral/atracacoes_paper.csv') %>% #local time
  mutate(start = dmy_hm(start),
         end = dmy_hm(end)) %>%
  subset(Natureza == "Granel Líquido" | 
           Natureza == "Granel Sólido" |
           Natureza == "Carga Geral" |
           Natureza == "Full Conteiner") %>%
  unique() #%>%
#select(Natureza, start, end)

atracacoes2 <- fread('./data/geral/atracacoes_paper2.csv') %>% #local time
  mutate(start = dmy_hm(start),
         end = dmy_hm(end)) %>%
  subset(Natureza == "Granel Líquido" | 
           Natureza == "Granel Sólido" |
           Natureza == "Carga Geral" |
           Natureza == "Full Conteiner") %>%
  unique() #%>%
# select(Natureza, start, end)

atracacoes <- bind_rows(atracacoes1, atracacoes2)

rm(atracacoes1, atracacoes2)


# Considerando Full Container == Carga Geral
atracacoes$Natureza[atracacoes$Natureza == "Full Conteiner"] <- "FULL CONTAINER"
atracacoes$Natureza[atracacoes$Natureza == "Carga Geral"] <- "FULL CONTAINER"
atracacoes$Natureza[atracacoes$Natureza == "Granel Sólido"] <- "SOLID BULK"
atracacoes$Natureza[atracacoes$Natureza == "Granel Líquido"] <- "LIQUID BULK"




#### 08/06/2025 --- 
duracao <- atracacoes %>%
  mutate(duracao = difftime(end, start, units = "hour"),
         duracao = round(as.numeric(gsub(" hours", "", duracao)), 0))


duracao <- atracacoes %>%
  mutate(
    duracao = as.numeric(difftime(end, start, units = "mins")), # duração em minutos
    duracao_15min = round(duracao / 15) # duração em blocos de 15 minutos
  )




#data wrangling: docking in/out by day

time <- seq(ymd_hms("2024-09-01 00:00:00"),
            ymd_hms("2024-12-04 23:45:00"),
            by="1 min") %>% as.data.frame()

colnames(time) <- c("date")

out <- vector("list", length(atracacoes))

for (i in 1:nrow(atracacoes)) {
  start <- min(atracacoes$start[i], atracacoes$end[i])
  end <- max(atracacoes$start[i], atracacoes$end[i])
  subperiodo <- subset(time, date <= (ceiling_date(end, unit = "1 min") + minutes(20)) & date >= (floor_date(start, unit = "1 min") - minutes(20)))
  out[[i]] <- subperiodo
  out[[i]]$Natureza <- atracacoes$Natureza[i]
  out[[i]]$Berco <- atracacoes$Berço[i]
  out[[i]]$start <- start
  out[[i]]$end <- end
  out[[i]]$moored <- ifelse((out[[i]]$date >= out[[i]]$start & out[[i]]$date <= out[[i]]$end), 1, 0)
  out[[i]]$un_berthing <- ifelse((
    (out[[i]]$date <= out[[i]]$start + minutes(20) &
       out[[i]]$date >= out[[i]]$start - minutes(20)) |
      (out[[i]]$date >= out[[i]]$end - minutes(20) &
         out[[i]]$date <= out[[i]]$end + minutes(20))), 1, 0)
}  


str(out)
output <- dplyr::bind_rows(out) %>% distinct() # combine the output into a single data frame



################################################################################
################################################################################
########## Data wrangling ---
################################################################################
################################################################################


# periodo 07/09/24 - 03/12/24

output <- output %>%
  mutate(#date = floor_date(date, "15 mins"),
         duracao_hrs = round(as.numeric(difftime(end, start, units = "hours")), 0)) 

#####################################################
################ DATASET #################

## Data Preparation ----
# 0. Baseline
# number of observations time (8448) * sensors(3) = 25344 * Natureza (3)= 76032

# Example lists
cargo_types <- c("FULL CONTAINER", "SOLID BULK", "LIQUID BULK")
sensor_height <- c("ALPHA", "BETA", "GAMMA")

baseline <- merge(sensor_height, cargo_types) |>
  merge(time) |>
  subset(date >= ymd_hms("2024-09-07 00:00:00") & 
           date <= ymd_hms("2024-12-03 23:45:00"))

colnames(baseline) <- c("site", "Natureza", "date")





# Suppose you already have a pollutants dataframe like:
# pollutants_df: with columns (cargo, sensor, pollutant, value)

# Join pollutants data to the grid (keeps structure)
df_15min <- baseline %>%
  left_join(pollutants, by = c("date", "site")) |>
  left_join(output, by = c('Natureza', 'date'))  %>%
  unique() |>
  mutate(date = floor_date(date, "15 mins")) |>
  group_by(site, Natureza, date) |>
  summarise(
    pm25 = mean(pm2p5, na.rm = T),
    so2 = mean(so2, na.rm = T),
    no2 = mean(no2, na.rm = T),
    moored = sum(moored, na.rm = T), # sum by minutes (i.e. minutes all vessels moored)
    n_vessels_time = sum(un_berthing, na.rm = T), # sum by minutes (i.e. minutes all vessels un_berthing)
    duracao_hrs = sum(duracao_hrs, na.rm = T))



df_20min <- baseline %>%
  left_join(pollutants, by = c("date", "site")) |>
  left_join(output, by = c('Natureza', 'date'))  %>%
  unique() |>
  mutate(date = floor_date(date, "20 mins")) |>
  group_by(site, Natureza, date) |>
  summarise(
    pm25 = mean(pm2p5, na.rm = T),
    so2 = mean(so2, na.rm = T),
    no2 = mean(no2, na.rm = T),
    moored = sum(moored, na.rm = T), # sum by minutes (i.e. minutes all vessels moored)
    n_vessels_time = sum(un_berthing, na.rm = T), # sum by minutes (i.e. minutes all vessels un_berthing)
    duracao_hrs = sum(duracao_hrs, na.rm = T))


df_1hour <- baseline %>%
  left_join(pollutants, by = c("date", "site")) |>
  left_join(output, by = c('Natureza', 'date'))  %>%
  unique() |>
  mutate(date = floor_date(date, "1 hour")) |>
  group_by(site, Natureza, date) |>
  summarise(
    pm25 = mean(pm2p5, na.rm = T),
    so2 = mean(so2, na.rm = T),
    no2 = mean(no2, na.rm = T),
    moored = sum(moored, na.rm = T), # sum by minutes (i.e. minutes all vessels moored)
    n_vessels_time = sum(un_berthing, na.rm = T), # sum by minutes (i.e. minutes all vessels un_berthing)
    duracao_hrs = sum(duracao_hrs, na.rm = T))





# ==================================================
# Separate by cargo datasets
# ==================================================


#1. Multiple Terminals  ---- 
# Assumption: do not matter the type of cargo, only



############################# 15 min ----

# Ensure factors
df_15min$site <- factor(df_15min$site)
df_15min$Natureza <- factor(df_15min$Natureza)




# Adding docking by cargo
df_multiples_bycargo_15min <- df_15min |>
  mutate(
    docking_day = as.factor(
      ifelse(n_vessels_time > 0,  ## make a field called "feature" and fill: make all values = "other"
      "TRUE", "FALSE")),
    docking_container = as.factor(
      ifelse(
        (n_vessels_time > 0 & Natureza == "FULL CONTAINER"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE")),
    docking_solids = as.factor(
      ifelse(
        (n_vessels_time > 0 & Natureza == "SOLID BULK"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE")),
    docking_liquids = as.factor(
      ifelse(
        (n_vessels_time > 0 & Natureza == "LIQUID BULK"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE"))
  ) 







############################# 20 min ----
# Ensure factors
df_20min$site <- factor(df_20min$site)
df_20min$Natureza <- factor(df_20min$Natureza)




# Adding docking by cargo
df_multiples_bycargo_20min <- df_20min |>
  mutate(
    docking_day = as.factor(
      ifelse(n_vessels_time > 0,  ## make a field called "feature" and fill: make all values = "other"
             "TRUE", "FALSE")),
    docking_container = as.factor(
      ifelse(
        (n_vessels_time > 0 & Natureza == "FULL CONTAINER"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE")),
    docking_solids = as.factor(
      ifelse(
        (n_vessels_time > 0 & Natureza == "SOLID BULK"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE")),
    docking_liquids = as.factor(
      ifelse(
        (n_vessels_time > 0 & Natureza == "LIQUID BULK"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE"))
  ) 









############################# 1 hour ----
# Ensure factors
df_1hour$site <- factor(df_1hour$site)
df_1hour$Natureza <- factor(df_1hour$Natureza)




# Adding docking by cargo
df_multiples_bycargo_1hour <- df_1hour |>
  mutate(
    docking_day = as.factor(
      ifelse(n_vessels_time > 0,  ## make a field called "feature" and fill: make all values = "other"
             "TRUE", "FALSE")),
    docking_container = as.factor(
      ifelse(
        (n_vessels_time > 0 & Natureza == "FULL CONTAINER"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE")),
    docking_solids = as.factor(
      ifelse(
        (n_vessels_time > 0 & Natureza == "SOLID BULK"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE")),
    docking_liquids = as.factor(
      ifelse(
        (n_vessels_time > 0 & Natureza == "LIQUID BULK"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE"))
  ) 



## Adding meteo ----
# WIND COMPONENTS: Create columns to and calculate the u and v wind components to all data ----
### a: ws / b: wd

u.wind <- function(a,b) {
  u <- - a * sin(2 * pi * b/360)
  return(u)
}


v.wind <- function(a,b) {
  v <- - a * cos(2 * pi * b/360)
  return(v)
}



# Calculate wind direction for this hour-long period (WD REAL)
wind_direction <- function(u,v) {
  x <- (atan2(u,v) * 360/2/pi) + 180
  return(x)
}


# Calculate vectorial wind speed for this hour-long period
wind_speed <- function(u,v) {
  z <- ((u^2 + v^2)^0.5)
  return(z)
}


# MODEL A - MANEUVERING EFFECT: 15 min ----

meteo_15min <- meteo  |>
  subset(date >= ymd_hms("2024-09-07 00:00:00") & 
           date <= ymd_hms("2024-12-03 23:45:00")) |>
  mutate(date = floor_date(date, "15 min"),
         u = u.wind(ws, wd),
         v = v.wind(ws, wd)) %>%
  group_by(date) |>
  summarise(
    ws = mean(ws, na.rm = TRUE),
    u.wind = mean(u, na.rm = TRUE),
    v.wind = mean(v, na.rm = TRUE)
  ) %>%
  mutate(wd = wind_direction(u.wind, v.wind)) %>%
  left_join(., meteo_simepar[, c("date", "temperatura", "umidade", "precipitacao")], by = 'date') 




#[ PRECISA CALCULAR MEDIAS AQUI E NO DE 20]
df_multiples_bycargo_15min <- df_multiples_bycargo_15min |>
  left_join(meteo_15min, by = "date")

colnames(df_multiples_bycargo_15min)




# MODEL A - MANEUVERING EFFECT: 20 min ----

meteo_20min <- meteo  |>
  full_join(time, by = c("date")) |>
  subset(date >= ymd_hms("2024-09-07 00:00:00") & 
           date <= ymd_hms("2024-12-03 23:45:00")) |>
  left_join(meteo_simepar[, c("date", "temperatura", "umidade", "precipitacao")], by = 'date') %>%
  mutate(date = floor_date(date, "20 min"),
         u = u.wind(ws, wd),
         v = v.wind(ws, wd)) %>%
  group_by(date) |>
  summarise(
    ws = mean(ws, na.rm = TRUE),
    u.wind = mean(u, na.rm = TRUE),
    v.wind = mean(v, na.rm = TRUE),
    temperatura = mean(temperatura, na.rm = TRUE),
    umidade = mean(umidade, na.rm = TRUE),
    precipitacao = mean(precipitacao, na.rm = TRUE)
  ) %>%
  mutate(wd = wind_direction(u.wind, v.wind))

df_multiples_bycargo_20min <- df_multiples_bycargo_20min |>
  left_join(meteo_20min, by = "date")

colnames(df_multiples_bycargo_20min)







# MODEL B - HOTELLING EFFECT ----

meteo_1hour <- meteo  |>
  subset(date >= ymd_hms("2024-09-07 00:00:00") & 
           date <= ymd_hms("2024-12-03 23:45:00")) |>
  mutate(date = floor_date(date, "1 hour"),
         u = u.wind(ws, wd),
         v = v.wind(ws, wd)) %>%
  left_join(., meteo_simepar[, c("date", "temperatura", "umidade", "precipitacao")], by = 'date') %>%
  group_by(date) |>
  summarise(
    ws = mean(ws, na.rm = TRUE),
    u.wind = mean(u, na.rm = TRUE),
    v.wind = mean(v, na.rm = TRUE),
    temperatura = mean(temperatura, na.rm = TRUE),
    umidade = mean(umidade, na.rm = TRUE),
    precipitacao = mean(precipitacao, na.rm = TRUE)
  ) %>%
  mutate(wd = wind_direction(u.wind, v.wind))



df_1hour$site <- factor(df_1hour$site)
df_1hour$Natureza <- factor(df_1hour$Natureza)



df_multiples_bycargo_1hour <- df_multiples_bycargo_1hour |>
  left_join(meteo_1hour, by = "date")  |>
  left_join(traj_10, by = "date")

colnames(df_multiples_bycargo_1hour)




# MODEL B - HOTELLING EFFECT : 1 day ----
df_multiples_bycargo_1day <- df_multiples_bycargo_15min  |>
  mutate(date = floor_date(date, "1 day")) %>%
  left_join(traj_10, by = "date") %>%
  group_by(site, Natureza, date) %>%
  summarise(
    ws = mean(ws, na.rm = TRUE),
    u = mean(u.wind, na.rm = TRUE),
    v = mean(v.wind, na.rm = TRUE),
    mixdepth = mean(mixdepth, na.rm = TRUE),
    temperatura = mean(temperatura, na.rm = TRUE),
    umidade = mean(umidade, na.rm = TRUE),
    precipitacao = sum(precipitacao, na.rm = TRUE),
    no2 = mean(no2, na.rm = TRUE),
    so2 = mean(so2, na.rm = TRUE),
    pm25 = mean(pm25, na.rm = TRUE),
    duracao_hrs = sum(duracao_hrs, na.rm = TRUE),
    n_vessels_time = sum(n_vessels_time, na.rm = TRUE)
  ) %>%
  mutate(wd = wind_direction(u, v),
         docking_container = as.factor(
           ifelse(
             (n_vessels_time > 0 & Natureza == "FULL CONTAINER"),  ## make a field called "feature" and fill: make all values = "other"
             "TRUE", "FALSE")),
         docking_solids = as.factor(
           ifelse(
             (n_vessels_time > 0 & Natureza == "SOLID BULK"),  ## make a field called "feature" and fill: make all values = "other"
             "TRUE", "FALSE")),
         docking_liquids = as.factor(
           ifelse(
             (n_vessels_time > 0 & Natureza == "LIQUID BULK"),  ## make a field called "feature" and fill: make all values = "other"
             "TRUE", "FALSE"))
  ) 



################################################################################
# Case-Crossover clogit Analysis: Advanced Sensitivity + Diagnostic Modules
################################################################################
# Features:
# - Sensitivity analysis for covariate sets (including meteorology and event duration)
# - Monte Carlo diagnostics (bias, variance, coverage, MCSE)
# - Strata diagnostics (usable strata, mean cases/controls per stratum)
# - Descriptions for interpreting diagnostics
################################################################################

library(dplyr)
library(lubridate)
library(survival)
library(broom)
library(readr)
library(ggplot2)
library(purrr)
library(openair)

##############################
### PARAMETERS AND SETUP   ###
##############################
pollutants    <- c("pm25", "no2", "so2")
cargo_types   <- c("SOLID BULK", "LIQUID BULK", "FULL CONTAINER")
sensor_types  <- c("ALPHA", "BETA", "GAMMA")
case_quantiles <- c(0.75, 0.90)
max_controls  <- 4
min_cases     <- 10
min_strata    <- 20

##############################
## Covariate Sensitivity Sets
##############################
# Model A: Berthing
cov_sets_A <- list(
  c("ws", "wd", "precipitacao"),                        # base
  c("ws", "wd", "precipitacao", "temperatura"),          # add temperatura
  c("ws", "wd", "precipitacao", "umidade"),              # add umidade
  c("ws", "wd", "precipitacao", "temperatura", "umidade"),# both
  c("ws", "wd", "precipitacao", "n_vessels_time"),                        # base
  c("ws", "wd", "precipitacao", "temperatura", "n_vessels_time"),          # add temperatura
  c("ws", "wd", "precipitacao", "umidade", "n_vessels_time"),              # add umidade
  c("ws", "wd", "precipitacao", "temperatura", "umidade", "n_vessels_time")# both
)

# Model B: Moored
cov_sets_B <- list(
  c("ws", "wd", "precipitacao", "duracao_hrs"),                               # base + duration
  c("ws", "wd", "precipitacao", "temperatura", "umidade", "duracao_hrs"),     # add temp/hum + duration
  c("ws", "wd", "precipitacao", "temperatura", "umidade", "mixdepth", "duracao_hrs") # all
)

##############################
### Helper Functions       ###
##############################

cap_controls <- function(dat, case_col = "case", stratum_col = "strata_id", max_controls = max_controls) {
  dat %>%
    group_by(.data[[stratum_col]]) %>%
    group_modify(~{
      cases <- .x %>% filter(.data[[case_col]] == 1)
      controls <- .x %>% filter(.data[[case_col]] == 0)
      controls_sampled <- controls %>% slice_sample(n = min(nrow(cases) * max_controls, nrow(controls)))
      bind_rows(cases, controls_sampled)
    }) %>%
    ungroup()
}

diagnostics_strata <- function(df_clogit, stratum_col = "strata_id") {
  strata_table <- df_clogit %>%
    group_by(.data[[stratum_col]]) %>%
    summarise(cases = sum(case == 1), controls = sum(case == 0), .groups = "drop")
  list(
    n_strata = n_distinct(df_clogit[[stratum_col]]),
    mean_cases_per_stratum = mean(strata_table$cases),
    mean_controls_per_stratum = mean(strata_table$controls),
    strata_table = strata_table
  )
}

montecarlo_diagnostics <- function(df, formula_str, n_sim = 50) {
  or_vec <- numeric(n_sim)
  for (i in 1:n_sim) {
    df_sim <- df %>% group_by(strata_id) %>% sample_frac(1, replace = TRUE) %>% ungroup()
    fit <- tryCatch(clogit(as.formula(formula_str), data = df_sim), error = function(e) NULL)
    if (!is.null(fit)) {
      summary_fit <- summary(fit)
      or_vec[i] <- summary_fit$coefficients[1,2] # log(OR)
    } else {
      or_vec[i] <- NA
    }
  }
  or_vec <- exp(or_vec)
  bias <- mean(or_vec, na.rm = TRUE) - 1
  variance <- var(or_vec, na.rm = TRUE)
  coverage <- mean(or_vec > 1, na.rm = TRUE)
  mcse <- sqrt(variance/n_sim)
  data.frame(Bias = bias, Variance = variance, Coverage = coverage, MCSE = mcse)
}

##############################
### Main Model Function     ###
##############################
# CHANGES IN MODEL B - case (no n_vessels_time criteria, since everytime there is one moored) and strata_id (no block, + month instead)

run_clogit_analysis <- function(
    df, model_type = c("A", "B"), 
    pollutants, cargo_types, sensor_types, case_quantile, covariate_set,
    max_controls, min_cases, min_strata
) {
  results.temp <- data.frame()
  diagnostics.temp <- list()
  
  for (poll in pollutants) {
    for (cargo in cargo_types) {
      for (sensor in sensor_types) {
        threshold <- quantile(df[[poll]][df$site == sensor], probs = case_quantile, na.rm = TRUE)
        if (is.na(threshold)) next
        
        # MODEL A: Berthing
        if (model_type == "A") {
          df_case <- df %>%
            mutate(hour = hour(date),
                   dow = wday(date),
                   month = month(date),
                   case = (.data[[poll]] > threshold) & Natureza == cargo & n_vessels_time > 0 & site == sensor,
                   strata_id = interaction(site, hour, dow, month, sep = "_", drop = TRUE))
          formula_str <- paste0("case ~ ", poll, " + ", paste(covariate_set, collapse = " + "), " + strata(strata_id)")
        }
        # MODEL B: Moored
        if (model_type == "B") {
          df_case <- df %>%
            mutate(hour = hour(date),
                   dow = wday(date),
                   month = month(date),
                   case = (.data[[poll]] > threshold) & Natureza == cargo & site == sensor,
                   #strata_week = floor_date(date, unit = "month"),
                   strata_id = interaction(site, hour, dow, month, sep = "_", drop = TRUE))
          formula_str <- paste0("case ~ ", poll, " + ", paste(covariate_set, collapse = " + "), " + strata(strata_id)")
        }
        
        df_case <- df_case %>% filter(!is.na(ws), !is.na(wd), !is.na(precipitacao))
        df_case <- df_case[order(df_case$strata_id),]
        
        df_clogit <- cap_controls(df_case, case_col = "case", stratum_col = "strata_id", max_controls = max_controls)
        df_clogit <- df_clogit %>%
          group_by(strata_id) %>%
          filter(any(case == 1), any(case == 0)) %>%
          ungroup()
        
        n_cases <- sum(df_clogit$case == 1, na.rm = TRUE)
        n_strata <- n_distinct(df_clogit$strata_id)
        n_covars <- length(covariate_set)
        if (n_cases < min_cases * n_covars | n_strata < min_strata) next
        
        diag <- diagnostics_strata(df_clogit)
        model_clogit <- tryCatch(clogit(as.formula(formula_str), data = df_clogit), error = function(e) NULL)
        if (is.null(model_clogit)) next
        
        summary_model <- summary(model_clogit)
        beta <- summary_model$coefficients[1,2]
        se <- summary_model$coefficients[1,3]
        beta_perc <- (beta - 1) * 100
        lo95ci <- (beta - 1.96 * se)
        hi95ci <- (beta + 1.96 * se)
        lo95ci_perc <- (lo95ci - 1) * 100
        hi95ci_perc <- (hi95ci - 1) * 100
        model_stats <- glance(model_clogit)
        
        mc_diag <- montecarlo_diagnostics(df_clogit, formula_str, n_sim = 50)
        
        results.temp <- bind_rows(results.temp, data.frame(
          Model = ifelse(model_type == "A", "Berthing", "Moored"),
          Covariates = paste(covariate_set, collapse = "+"),
          Strata = paste0(ifelse(model_type == "A", "site_hour_dow_month", "site_hour_dow_month")),
          Case = paste0("Qtil", case_quantile * 100, "_sensor"),
          Outcome = paste0(toupper(poll), " Episode"),
          Group = cargo,
          Elevation = sensor,
          Source = "Shipping + Weather",
          Controls_n = sum(df_clogit$case == 0, na.rm = TRUE),
          Cases_n = sum(df_clogit$case == 1, na.rm = TRUE),
          OddsRatio = beta,
          Lower_95CI = lo95ci,
          Upper_95CI = hi95ci,
          Perc_increase = beta_perc,
          Perc_increase_Lower_95CI = lo95ci_perc,
          Perc_increase_Upper_95CI = hi95ci_perc,
          AIC = model_stats$AIC,
          BIC = model_stats$BIC,
          poll_pvalue = summary_model$coefficients[1, "Pr(>|z|)"],
          MC_Bias = mc_diag$Bias,
          MC_Variance = mc_diag$Variance,
          MC_Coverage = mc_diag$Coverage,
          MC_MCSE = mc_diag$MCSE,
          Strata_n = diag$n_strata,
          MeanCases_perStrata = diag$mean_cases_per_stratum,
          MeanControls_perStrata = diag$mean_controls_per_stratum
        ))
        diagnostics.temp[[paste(poll, cargo, sensor, paste(covariate_set, collapse = "+"), sep = "_")]] <- list(
          strata_diag = diag,
          montecarlo_diag = mc_diag,
          formula = formula_str
        )
      }
    }
  }
  return(list(results = results.temp, diagnostics = diagnostics.temp))
}

#! CHANGES IN MODEL B - case (no n_vessels_time criteria, since everytime there is one moored) and strata_id (no block, + month instead)


##############################
### Sensitivity Analysis    ###
##############################
run_sensitivity <- function(df, model_type, pollutants, cargo_types, sensor_types,
                            quantiles, cov_sets, max_controls, min_cases, min_strata) {
  all_results <- data.frame()
  all_diags <- list()
  for (case_quantile in quantiles) {
    for (covariate_set in cov_sets) {
      out <- run_clogit_analysis(
        df = df,
        model_type = model_type,
        pollutants = pollutants,
        cargo_types = cargo_types,
        sensor_types = sensor_types,
        case_quantile = case_quantile,
        covariate_set = covariate_set,
        max_controls = max_controls,
        min_cases = min_cases,
        min_strata = min_strata
      )
      all_results <- bind_rows(all_results, out$results)
      all_diags[[paste0("Q", case_quantile*100, "_", paste(covariate_set, collapse = "+"))]] <- out$diagnostics
    }
  }
  return(list(results = all_results, diagnostics = all_diags))
}

##############################
### RUN FULL ANALYSIS       ###
##############################



#0. [NOT GOOD][NOT GOOD] modelA testing SO2 for hourly 1day, no hour strata, qtil75, qtil90
#1. [NO DIFFERENCE FOR SO2 GAMMA] modelA testing SO2 for hourly 15 min, hour strata, qtil50 ## all_results <- results_A$results |> subset(Case == "Qtil50_sensor") 
#2. [NOT GOOD]modelA testing SO2 for hourly 15 min, no hour strata
#3. [ ok] testing adding n_vessels_time to berth model

results_A <- run_sensitivity(
  #df = df_multiples_bycargo_15min,
  df = df_multiples_bycargo_20min, #primary, corrected case by Natureza too
  #df = df_multiples_bycargo_1hour, #sensitivity
  model_type = "A",
  pollutants = pollutants,
  cargo_types = cargo_types,
  sensor_types = sensor_types,
  quantiles = case_quantiles,
  cov_sets = cov_sets_A,
  max_controls = max_controls,
  min_cases = min_cases,
  min_strata = min_strata
)
write_csv(results_A$results, "./output/results_case_crossover_A20min_sens.csv")

results_B <- run_sensitivity(
  df = df_multiples_bycargo_1hour,
  #df = df_multiples_bycargo_20min,
  model_type = "B",
  pollutants = pollutants,
  cargo_types = cargo_types,
  sensor_types = sensor_types,
  quantiles = case_quantiles,
  cov_sets = cov_sets_B,
  max_controls = max_controls,
  min_cases = min_cases,
  min_strata = min_strata
)
write_csv(results_B$results, "./output/results_case_crossover_B1hour_sens.csv")


# Suppose results_A and results_B are outputs from run_sensitivity (each is a list with $results)
# Each dataframe already has Model and Covariates columns from the code above.
results_A <- read_csv("./output/results_case_crossover_A20min_sens.csv")
results_B <- read_csv("./output/results_case_crossover_B1hour_sens.csv")
all_results <- bind_rows(results_A$results, results_B$results)

# Bind results together ----
all_results <- bind_rows(results_A$results, results_B$results)

View(results_A$results)


# If you want to ensure/standardize the identifiers:
all_results <- all_results %>%
  mutate(
    ModelType = Model,
    CovariateSet = Covariates
  )

# Save as CSV
write_csv(all_results, "./output/results_casecrossover_all_full_primary.csv")

# Optionally, view or explore
View(all_results)

all_results <- fread("./output/results_casecrossover_all_full_primary.csv")


##############################
### Diagnostics Guide      ###
##############################
diagnostics_description <- "
Diagnostics parameters:
- OddsRatio: Estimated odds ratio for pollutant effect; >1 means higher odds, <1 lower odds.
- Lower_95CI, Upper_95CI: 95% confidence interval for OddsRatio.
- Perc_increase: Percent change in odds per unit pollutant increase.
- AIC, BIC: Model fit statistics; lower is better.
- MC_Bias: Monte Carlo estimate of bias (mean deviation from null, i.e., 1.0 for OR).
- MC_Variance: Monte Carlo variance of odds ratio.
- MC_Coverage: Proportion of MC simulations where OR > 1 (if true effect is >1, higher is better).
- MC_MCSE: Monte Carlo standard error (precision of MC estimate).
- Strata_n: Number of usable strata in the model (should be >= min_strata).
- MeanCases_perStrata, MeanControls_perStrata: Average cases/controls per stratum (higher is better for stability).
Interpretation:
- Prefer models with higher case/control counts, more usable strata, lower AIC/BIC, and stable Monte Carlo diagnostics (low bias, MCSE).
- If MC_Coverage is near 0.95, confidence intervals are well-calibrated.
- If MC_Bias is near 0, model is unbiased.
"

cat(diagnostics_description)


#######################{ BEST DIAGNOSTIC BY "BLOCK" }#######################

library(dplyr)
library(ggplot2)
library(readr)
library(dplyr)
library(ggplot2)
library(readr)

# 1. Define block and spec columns
block_cols <- c("Outcome", "Group", "Elevation")
spec_cols <- c("Covariates", "Strata", "Case")

# 2. Count blocks per model specification
valid_results <- all_results %>%
  filter(Cases_n >= 10, Controls_n >= 10, Strata_n >= 20,
         !is.na(AIC), !is.na(MC_Bias), !is.na(MC_MCSE)) %>%
  mutate(abs_MC_Bias = abs(MC_Bias),
         abs_MC_MCSE = abs(MC_MCSE),
         Coverage_dist = abs(MC_Coverage - 0.95),
         is_signif = poll_pvalue < 0.05)

block_counts <- valid_results %>%
  group_by(across(all_of(spec_cols))) %>%
  summarise(
    n_blocks = n(),   # number of Outcome × Group × Elevation combos with valid results
    mean_AIC = mean(AIC, na.rm=TRUE),
    mean_MC_Bias = mean(abs_MC_Bias, na.rm=TRUE),
    mean_MC_MCSE = mean(abs_MC_MCSE, na.rm=TRUE),
    mean_Coverage_dist = mean(Coverage_dist, na.rm=TRUE),
    n_signif = sum(is_signif, na.rm=TRUE),
    .groups = "drop"
  )

# 3. Select the spec with the highest n_blocks, then best diagnostics
best_global_spec <- block_counts %>%
  arrange(desc(n_blocks), mean_AIC, mean_MC_Bias, mean_MC_MCSE, mean_Coverage_dist, desc(n_signif)) %>%
  slice(1)

# 4. Filter all valid results for this best specification (for both Model A and B)
plot_best <- valid_results %>%
  filter(Covariates == best_global_spec$Covariates,
         Strata == best_global_spec$Strata,
         Case == best_global_spec$Case)

# 5. Plot for Model A (Berthing)
# Okabe-Ito palette (3 colors, colorblind safe)
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73") # orange, sky blue, bluish green

plot_A_best <- plot_best %>% filter(Model == "Berthing")
pa <- ggplot(
  plot_A_best %>%
    mutate(
      Elevation = factor(Elevation, levels = c("ALPHA", "BETA", "GAMMA"),
                         labels = c("ALPHA \n (1.5 m AGL)", "BETA \n (10 m AGL)", "GAMMA \n (40 m AGL)")),
      Outcome = factor(Outcome, levels = c("SO2 Episode", "NO2 Episode", "PM25 Episode"),
                       labels = c(expression(SO[2]~"Episode"), expression(NO[2]~"Episode"), expression(PM[2.5]~"Episode"))),
      Group = factor(Group, levels = c("LIQUID BULK", "SOLID BULK", "FULL CONTAINER"),
                     labels = c("LIQUID BULK", "SOLID BULK", "FULL CONTAINER"))
    ),
  aes(y = Perc_increase, x = Elevation,
      ymin = Perc_increase_Lower_95CI,
      ymax = Perc_increase_Upper_95CI,
      color = Group)
) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = okabe_ito) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "Elevation", y = "Percentage change (95% CI)",
       title = paste("Best Global Model Spec (Berthing):", best_global_spec$Covariates, "|", best_global_spec$Strata, "|", best_global_spec$Case)) +
  theme_classic() +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(. ~ Outcome, labeller = label_parsed, scales = "free")



# 7. Save results for reporting
write_csv(plot_best, "./output/results_casecrossover_best_global_modelA.csv")

################################################################################
# This script finds the unique best model specification (Covariates, Strata, Case)
# that gives valid results for the most blocks (pollutant × cargo × sensor),
# and among those, selects the best diagnostics.
# It then plots all valid results for that model specification for Model A.
################################################################################



################################################################################
# Select Best Model B Spec (Moored) Regardless of Model A, by Block Coverage and Diagnostics
################################################################################
# Finds the best Model B specification (Covariates, Strata, Case) that yields
# valid results for the greatest number of pollutant × cargo × sensor blocks,
# then selects among those the best diagnostics (AIC, MC_Bias, etc.).
# Plots all valid results for that specification. Robust to empty data.
################################################################################

library(dplyr)
library(ggplot2)
library(readr)

# 1. Define block and spec columns
block_cols <- c("Outcome", "Group", "Elevation")
spec_cols <- c("Covariates", "Strata", "Case")

# 2. Filter only Model B (Moored) valid results
valid_results_B <- all_results %>%
  filter(Model == "Moored",
         Cases_n >= 10, Controls_n >= 10, Strata_n >= 20,
         !is.na(AIC), !is.na(MC_Bias), !is.na(MC_MCSE)) %>%
  mutate(abs_MC_Bias = abs(MC_Bias),
         abs_MC_MCSE = abs(MC_MCSE),
         Coverage_dist = abs(MC_Coverage - 0.95),
         is_signif = poll_pvalue < 0.05)

# 3. For each Model B spec, count number of blocks it covers and summarize diagnostics
block_counts_B <- valid_results_B %>%
  group_by(across(all_of(spec_cols))) %>%
  summarise(
    n_blocks = n(),   # number of Outcome × Group × Elevation combos with valid results
    mean_AIC = mean(AIC, na.rm=TRUE),
    mean_MC_Bias = mean(abs_MC_Bias, na.rm=TRUE),
    mean_MC_MCSE = mean(abs_MC_MCSE, na.rm=TRUE),
    mean_Coverage_dist = mean(Coverage_dist, na.rm=TRUE),
    n_signif = sum(is_signif, na.rm=TRUE),
    .groups = "drop"
  )

# 4. Select the spec with the highest n_blocks, then best diagnostics
best_global_spec_B <- block_counts_B %>%
  arrange(desc(n_blocks), mean_AIC, mean_MC_Bias, mean_MC_MCSE, mean_Coverage_dist, desc(n_signif)) %>%
  slice(1)

# 5. Filter all valid Model B results for this best specification
plot_B_best <- valid_results_B %>%
  filter(Covariates == best_global_spec_B$Covariates,
         Strata == best_global_spec_B$Strata,
         Case == best_global_spec_B$Case)



# 6. Plot for Model B (Moored) - robust to empty data
if (nrow(plot_B_best) == 0) {
  message("No valid Model B results for the globally best specification.")
} else {
 pb <- ggplot(
    plot_B_best %>%
      mutate(
        Elevation = factor(Elevation, levels = c("GAMMA", "BETA", "ALPHA"),
                           labels = c("'GAMMA \n (40 m AGL)'", "'BETA \n (10 m AGL)'", "'ALPHA \n (1.5 m AGL)'")),
        Outcome = factor(Outcome, levels = c("SO2 Episode", "NO2 Episode", "PM25 Episode"),
                         labels = c(expression(SO[2]~"Episode"), expression(NO[2]~"Episode"), expression(PM[2.5]~"Episode"))),
        Group = factor(Group, levels = c("LIQUID BULK", "SOLID BULK", "FULL CONTAINER"),
                       labels = c("LIQUID \n BULK", "SOLID \n BULK", "FULL \n CONTAINER"))
      ),
    aes(y = Perc_increase, x = Group,
        ymin = Perc_increase_Lower_95CI,
        ymax = Perc_increase_Upper_95CI,
        color = Group)
  ) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = okabe_ito) +
    geom_errorbar(width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(x = "Cargo type", y = "Percentage change (95% CI)",
         title = paste("Best Global Model B Spec (Moored):", best_global_spec_B$Covariates, "|", best_global_spec_B$Strata, "|", best_global_spec_B$Case)) +
    theme_classic() +
    theme(text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(Elevation ~ Outcome, labeller = label_parsed)
}



# 7. Save results for reporting (if not empty)
if (nrow(plot_B_best) > 0) {
  write_csv(plot_B_best, "./output/results_casecrossover_best_global_modelB.csv")
}




# > Option 2 < #



pb <- ggplot(
  plot_B_best %>%
    mutate(
      Elevation = factor(Elevation, levels = c("ALPHA", "BETA", "GAMMA"),
                         labels = c("ALPHA \n (1.5 m AGL)", "BETA \n (10 m AGL)", "GAMMA \n (40 m AGL)")),
      Outcome = factor(Outcome, levels = c("SO2 Episode", "NO2 Episode", "PM25 Episode"),
                       labels = c(expression(SO[2]~"Episode"), expression(NO[2]~"Episode"), expression(PM[2.5]~"Episode"))),
      Group = factor(Group, levels = c("LIQUID BULK", "SOLID BULK", "FULL CONTAINER"),
                     labels = c("LIQUID BULK", "SOLID BULK", "FULL CONTAINER"))
    ),
  aes(y = Perc_increase, x = Elevation,
      ymin = Perc_increase_Lower_95CI,
      ymax = Perc_increase_Upper_95CI,
      color = Group)
) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = okabe_ito) +
  geom_errorbar(width = 0.2, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "Cargo type", y = "Percentage change (95% CI)",
       title = paste("Best Global Model B Spec (Moored):", best_global_spec_B$Covariates, "|", best_global_spec_B$Strata, "|", best_global_spec_B$Case)) +
  theme_classic() +
  theme(text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(. ~ Outcome, labeller = label_parsed, scales = "free")




################################################################################
# This script finds the unique best Model B specification (Covariates, Strata, Case)
# that gives valid results for the most blocks (pollutant × cargo × sensor) for Model B,
# and among those, selects the best diagnostics.
# It then plots all valid results for that model specification for Model B.
################################################################################





# Save as svg
library(patchwork)
library(Cairo)


svg(file = paste0("./output/figs/casecrossover_bestmodels_A20min_B1hour_facetwrap.svg"),
    width = 11, height = 8)



(pa / pb) + plot_layout(widths = c(1), heights = c(1, 1))


dev.off()










# END SCRIPT # -----






