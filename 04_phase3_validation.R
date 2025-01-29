library(pacman)
p_load("devtools", "scales", "ellipse", "lazyeval", "igraph",  "ggraph",
       "reshape2", "knitr", "stringr", "jsonlite", "rstudioapi", "tidyverse",
       "dampack", "data.table", "tornado", "ggplot2", "viridis",
       "kableExtra")


## Load model parameters
source("01_phase3_inputs.R")
l_params_all <- load_phase3_params()
v_names_str <- c("Status Quo", "Testing", "Vaccination", "Vaccination and Testing")
n_str <- length(v_names_str)

## Load model
source("02_phase3_functions.R")

## Run model
df_results  <- parvo_model(l_params_all)


