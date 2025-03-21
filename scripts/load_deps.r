## Load them all

library("tidymodels")
library("tidyverse")
library("embed")
library("feather")
library("stringr")
library("foreach")
library("glue")
library("doParallel")
library("future")
library("embed")
library("readxl")
library("biomaRt")

library("dotenv")
load_dot_env()

## get some env vars
proj_folder <- Sys.getenv("PROJ_FOLDER")
data_folder <- Sys.getenv("DATA_FOLDER")
raw_data_folder <- Sys.getenv("RAW_DATA_FOLDER")

## register BiocParallel cores
if(.Platform$OS.type == "windows") {
  makeCluster(4) %>% doParallel::registerDoParallel(4)
} else {
  doParallel::registerDoParallel(cores = 4)
}

## load the commonly used backend functions
source("./scripts/backend_env_functions.r")