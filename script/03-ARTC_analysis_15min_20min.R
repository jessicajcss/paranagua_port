################################################################################
# Script Name: Port Docking and Air Pollution Analysis (ART-C Statistical Tests)
#
# Purpose:
# This script processes air quality observations (PM2.5, SO2, NO2) and port
# vessel docking activity (by cargo type and terminal). It aggregates pollutant
# concentrations to 20-minute intervals, classifies docking events, and applies
# Aligned Rank Transform ANOVA (ART) followed by Tukey-adjusted pairwise
# contrasts to evaluate differences in pollutant levels:
#   (1) Between monitoring sites
#   (2) Between docking vs. non-docking periods
#   (3) Within combinations of site × cargo type.
#
# Inputs:
#   - ./data/geral/pollutants.csv           : Air pollutant time series data
#   - ./data/geral/atracacoes_paper.csv     : Vessel docking log (part 1)
#   - ./data/geral/atracacoes_paper2.csv    : Vessel docking log (part 2)
#
# Outputs:
#   - ./output/stats/pairwise_ARTC_resultsFOR*.csv : Pairwise contrasts per cargo group
#   - ./output/stats/pairwise_all_results_ARTC_20min.csv : Combined results table
#
# How to Run:
#   1. Ensure the required R packages are installed (see libraries below).
#   2. Place input files in the paths listed above.
#   3. Adjust analysis period in the `subset(date >= ... date <= ...)` sections if needed.
#   4. Run the script from start to finish in R or RStudio.
#
# Author: Santos-Silva, J. C.
# Last Updated: 05/09/2025
################################################################################


library(lubridate)
library(bazar)
library(rstatix)
library(dplyr)
library(tidyverse)
library(openair)



################################################################################
################################################################################
########## Previous run: Data INPUT ---
################################################################################
################################################################################


#==== POLLUTANTS ----
pollutants <- fread('./data/geral/pollutants.csv') %>% # time TZ; pm in ug/m³ & gases in ppb
  mutate(date = dmy_hm(date),
         site = case_when(sensor == "alfa" ~ "ALPHA", 
                          sensor == "gamma" ~ "BETA",
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
    duracao_15min = round(duracao / 20) # duração em blocos de 20 minutos
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



summary(df)

# ==================================================
# Separate by cargo datasets
# ==================================================

# 0. Select maneuver time interval = *_15min or *_20min  window

df <- df_15min # if 15 min window, be sure to correct the [out] list
df <- df_20min


# Ensure factors
df$site <- factor(df$site)
df$Natureza <- factor(df$Natureza)


#1. Multiple Terminals ----
library(dplyr)
library(tidyr)


df_multiples <- df  |>
  group_by(site, Natureza, date) |>
  summarise(
    pm25 = mean(pm25, na.rm = T),
    so2 = mean(so2, na.rm = T),
    no2 = mean(no2, na.rm = T),
    n_vessels_day = sum(n_vessels_time, na.rm = T),
    duracao_hrs = sum(duracao_hrs, na.rm = T)
  ) %>%
  mutate(docking_day = as.factor(ifelse(n_vessels_day > 0,  ## make a field called "feature" and fill: make all values = "other"
                                        "TRUE", "FALSE"))) 


# Ensure factors
df_multiples$docking_day <- factor(df_multiples$docking_day)



# Adding docking by cargo ----


df_container  <-  df_multiples |>
  mutate(
    docking_day = as.factor(
      ifelse(
        (n_vessels_day > 0 & Natureza == "FULL CONTAINER"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE")))

df_solids  <-  df_multiples |>
  mutate(
    docking_day = as.factor(
      ifelse(
        (n_vessels_day > 0 & Natureza == "SOLID BULK"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE")))

df_liquids  <-  df_multiples |>
  mutate(
    docking_day = as.factor(
      ifelse(
        (n_vessels_day > 0 & Natureza == "LIQUID BULK"),  ## make a field called "feature" and fill: make all values = "other"
        "TRUE", "FALSE")))




df <- df_multiples





################################################################################
################################################################################
########## ART Analysis ---
################################################################################
################################################################################
# ---------------------------------------
# Load necessary libraries
# ---------------------------------------
library(tidyverse)
library(rstatix)
library(ggpubr)
library(emmeans)
library(ARTool)
library(patchwork)
library(Cairo)


###########################################
################ TESTE 1 ##################
######## comparing docking among sites ----

# ==================================================

# ==================================================
# Function: Run ART-C (main and interaction contrasts) for a dataset, using Tukey p-value adjustment
# ==================================================
library(tidyverse)
library(ARTool)
library(emmeans)
library(readr)

run_artc_for_dataset <- function(df, cargo_name, pollutants, docking_var = "docking_day") {
  results_all <- list()
  
  for (pollutant in pollutants) {
    cat("\n---", cargo_name, "|", pollutant, "---\n")
    
    df_sub <- df %>%
      select(site, docking = !!sym(docking_var), value = all_of(pollutant)) %>%
      drop_na()
    
    # Fit ART model once
    art_model <- art(value ~ site * docking, data = df_sub)
    
    # -------------------------------------------------
    # MAIN EFFECT: Site (averaged over docking), using artlm.con() and Tukey adjustment
    # -------------------------------------------------
    emm_site <- emmeans(artlm.con(art_model, "site"), ~ site | docking)
    contrast_site <- contrast(emm_site, method = "pairwise", adjust = "tukey") %>%
      as_tibble() %>%
      mutate(
        pollutant = pollutant,
        test = "Site main effect (averaged over docking)",
        cargo_docking = cargo_name,
        p.adj = p.value # Tukey-adjusted p-value is returned by emmeans
      )
    
    # -------------------------------------------------
    # MAIN EFFECT: Docking (averaged over site), using artlm.con() and Tukey adjustment
    # -------------------------------------------------
    emm_docking <- emmeans(artlm.con(art_model, "docking"), ~ docking | site)
    contrast_docking <- contrast(emm_docking, method = "pairwise", adjust = "tukey") %>%
      as_tibble() %>%
      mutate(
        pollutant = pollutant,
        test = "Docking main effect (averaged over site)",
        cargo_docking = cargo_name,
        p.adj = p.value
      )
    
    # -------------------------------------------------
    # INTERACTION CONTRASTS: Site differences within Docking (ART-C) using art.con() with Tukey adjustment
    # -------------------------------------------------
    artc_site_docking <- art.con(art_model, "site:docking", adjust = "tukey") %>%
      as_tibble() %>%
      separate(contrast, into = c("group1","group2"), sep = " - ") %>%
      separate(group1, into = c("site1","docking1"), sep = ",") %>%
      separate(group2, into = c("site2","docking2"), sep = ",") %>%
      mutate(
        docking = docking1,
        group1 = site1,
        group2 = site2,
        pollutant = pollutant,
        test = "Site differences within docking (ART-C)",
        cargo_docking = cargo_name,
        p.adj = p.value
      )
    
    # -------------------------------------------------
    # INTERACTION CONTRASTS: Docking differences within Site (ART-C), using art.con() on swapped model with Tukey adjustment
    # -------------------------------------------------
    art_model_swapped <- art(value ~ docking * site, data = df_sub)
    artc_docking_site <- art.con(art_model_swapped, "docking:site", adjust = "tukey") %>%
      as_tibble() %>%
      separate(contrast, into = c("group1","group2"), sep = " - ") %>%
      separate(group1, into = c("docking1","site1"), sep = ",") %>%
      separate(group2, into = c("docking2","site2"), sep = ",") %>%
      mutate(
        site = site1,
        group1 = docking1,
        group2 = docking2,
        pollutant = pollutant,
        test = "Docking differences within site (ART-C)",
        cargo_docking = cargo_name,
        p.adj = p.value
      )
    
    # Combine all results for this pollutant
    results_all[[pollutant]] <- bind_rows(
      contrast_site,
      contrast_docking,
      artc_site_docking,
      artc_docking_site
    )
  }
  
  bind_rows(results_all)
}

# ==================================================
# Run for all cargo datasets
# ==================================================
cargo_list <- list(
  "Multiple Terminals"   = df,
  "Solid Bulk Terminal"  = df_solids,
  "Liquid Bulk Terminal" = df_liquids,
  "Container Terminal"   = df_container
)


pollutants <- c("so2", "no2", "pm25")
output_dir <- "./output/stats"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cargo_map <- list(
  "Multiple Terminals"   = "port",
  "Solid Bulk Terminal"  = "solids",
  "Liquid Bulk Terminal" = "liquids",
  "Container Terminal"   = "container"
)

all_results_list <- list()


for (cargo_name in names(cargo_list)) {
  res <- run_artc_for_dataset(cargo_list[[cargo_name]], cargo_name, pollutants)
  suffix <- cargo_map[[cargo_name]]
  out_file <- file.path(output_dir, paste0("pairwise_ARTC_resultsFOR", suffix, ".csv"))
  write_csv(res, out_file)
  all_results_list[[cargo_name]] <- res
}

# ==================================================
# Combine all into one big table
# ==================================================
all_results <- bind_rows(all_results_list, .id = "cargo_list_name")
write_csv(all_results, file.path(output_dir, "pairwise_all_results_ARTC_20min.csv"))

message("✅ Done. ART-C main and interaction contrasts (Tukey-adjusted p-values) written to: ", normalizePath(output_dir))
