#############################################################################
#############################################################################
# # Script para trajectory analysis
#   Production of HYSPLIT trajectory files
### Date: 21-01-2022
### Last update: 12-07-2025
### Elaborado por: Jéssica C. dos Santos-Silva

# # Publicação associada: xxx
#############################################################################
#############################################################################



## Referencias ---- 
#
# Production of HYSPLIT trajectory files:
#https://bookdown.org/david_carslaw/openair/sec-prod-hyspl-traj.html

# Trajectory analysis
#https://bookdown.org/david_carslaw/openair/sec-trajPlot.html


##############################################
############### LOAD PACKAGES ################

library(lubridate) # gerenciar dados data/hora 
library(tidyverse)
library(openair)
library(lattice)
library(devtools)
library(devEMF)

# aumentar memoria disponível para rodar dados
memory.size() 
memory.limit()

memory.limit(size = 28000)




## OBSERVAÇÕES ----

# https://openair-project.github.io/book/sections/appendices/appendix-hysplit.html

#1. Download and install the NOAA Hysplit model, somewhere with write access (see below).
#2. Download the monthly meteorological (.gbl) files also from the NOAA website.
#3. Obtain the code to run Hysplit.


##############################################
############### IMPORTAR DADOS ###############
##############################################


########
# Production of HYSPLIT trajectory files ----
########



#1. Download and install the NOAA Hysplit model, somewhere with write access


#2. Download the monthly meteorological (.gbl) files also from the NOAA website.
## corrected - renomear com RP+'name' depois

# Option 1 - SALVA NO DIRETORIO PRINCIPAL
getMet <- function (year = 2022:2025, month = 1:12) {
  for (i in seq_along(year)) {
    for (j in seq_along(month)) {
      download.file(url = paste0("ftp://arlftp.arlhq.noaa.gov/pub/archives/reanalysis/RP",
                                 year[i], sprintf("%02d", month[j]), ".gbl"),
                    destfile = paste0("RP", year[i],
                                      sprintf("%02d", month[j]), ".gbl"), mode = "wb")
    }
  }
}


###======
# Option 2 - TESTAR
getMet <- function(year, 
                   month, 
                   dest_folder) {
  # Create destination folder if it doesn't exist
  if (!dir.exists(dest_folder)) {
    dir.create(dest_folder, recursive = TRUE)
  }
  
  for (i in seq_along(year)) {
    for (j in seq_along(month)) {
      file_name <- paste0("RP", year[i], sprintf("%02d", month[j]), ".gbl")
      destfile <- file.path(dest_folder, file_name)
      
      # Skip download if file already exists
      if (!file.exists(destfile)) {
        url <- paste0("ftp://arlftp.arlhq.noaa.gov/pub/archives/reanalysis/", file_name)
        tryCatch({
          download.file(url = url, destfile = destfile, mode = "wb")
          message("Downloaded: ", file_name)
        }, error = function(e) {
          warning("Failed to download ", file_name, ": ", e$message)
        })
      } else {
        message("File already exists, skipping: ", file_name)
      }
    }
  }
}




####========


options(timeout = 300)
# Option 1
getMet(year = 2022, month = 1:12)
getMet(year = 2024, month = 8:12)

# Option 2
getMet(year = 2020,
       month = 1:12, 
       dest_folder = "./output/Hysplit/TrajData")


#3. Obtain the code to run Hysplit.


source_gist("https://gist.github.com/davidcarslaw/c67e33a04ff6e1be0cd7357796e4bdf5",
            filename = "run_hysplit.R")



# MODIFY THE function: parse_date_arguments with the command 
trace("parse_date_arguments", edit=TRUE)

# Edit date <- lubridate::parse_date_time(date, c("ymd", "dmy", "dmy_hm", "ymd_hm")) to
# date <- lubridate::parse_date_time(date, c("ymd", "dmy", "dmy_HM", "ymd_HM"))



# baixar trajetórias para os pontos amostrais
# SEE !!!it gives an error regarding column V1 read)
trajetorias_set <- trajetorias # check if the results given by the new function (disregarding v1) matches the original

source("./scripts/run_hysplit_sourcefile.R")
trajetorias <- run_hysplit( #713 TRAJECTORIES
  latitude = -25.50335, 
  longitude =  -48.5191, 
  runtime = -96, 
  start_height = 10, 
  model_height = 10000, 
  interval = "1 hour", 
  start = "2024-09-07",
  end = "2024-12-04",
  hysplit_exec = "C:/hysplit/exec", 
  hysplit_input = here::here("output/Hysplit/TrajData"), 
  hysplit_output = here::here("output/Hysplit/hysplit_output"),
  site = "paranagua")

#! trajetorias_set <- trajetorias # check if the results given by the new function (disregarding v1) matches the original

saveRDS(trajetorias, 
       #  file = here::here("output/Hysplit/hysplit_output/paranagua_forward48h_10m_1h.rds"))
       file = here::here("output/Hysplit/hysplit_output/paranagua_forward48h_40m_1h.rds"))
       # file = here::here("output/Hysplit/hysplit_output/paranagua_96h_10m_3h.rds"))
       # file = here::here("output/Hysplit/hysplit_output/paranagua_96h_100m_3h.rds"))
       #  file = here::here("output/Hysplit/hysplit_output/paranagua_96h_300m_3h.rds"))
       # file = here::here("output/Hysplit/hysplit_output/paranagua_96h_40m_3h.rds"))

# verificar número de horas de cada trajetoria
contagem <- table(as.Date(trajetorias$date)) #OK N=1176
View(contagem)

