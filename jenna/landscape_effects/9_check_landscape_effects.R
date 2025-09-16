##### Test accuracy of landscape effects #####
# Script Initiated: June 15, 2025
# By: Jenna Melanson
# Goals:
### Check whether models can accurately return estimates of landscape effects on foraging distance;
### assuming that a "true" landscape metric value is known for each 


##### Load packages #####
library(rstan)
library(matrixStats)
library(sp)
library(gstat)
library(ggplot2)
library(reshape2)
library(raster)
library(rasterVis)
library(parallel)
library(future)
library(furrr)
library(dplyr)
library(tidyr)
library(tibble)
library(posterior)
library(terra)

##### Set Environment #####
setwd("/Users/jenna1/Documents/UBC/bayes2025")

##### Source functions #####
source("jenna/src/GeneralizedSimFunctions.R")

##### Load in floral quality landscape #####
fq = readRDS("jenna/sim_data/random_field_range10/landscape_001.rds")
fq = terra::rast(fq)

##### Load in landscape metrics #####
config = readRDS("jenna/sim_data/random_field_range400/landscape_001.rds")
config = terra::rast(config)

##### Simulate some IDEAL data with landscape effects on foraging distance #####
result_ideal = draw_bees_colony_restricted(
  sample_size     = 1000,
  landscape_size  = 1500,
  colonygrid_size = 700,
  trapgrid_size   = 300,
  number_traps    = 25,
  number_colonies = 200,
  colony_sizes    = rep(100, 200),
  rho            = 50,
  theta           = 0.2,
  alpha = 0.5,
  configuration = config,
  resource_landscape = fq,
  nesting_landscape = NULL,
  distance_decay = "exponential"
)


##### Simulate some REALISTIC data with landscape effects on foraging distance #####
result_real = draw_bees_colony_restricted(
  sample_size     = 1000,
  landscape_size  = 1500,
  colonygrid_size = 700,
  trapgrid_size   = 300,
  number_traps    = 25,
  number_colonies = 2000,
  colony_sizes    = rep(100, 2000),
  rho            = 50,
  theta           = 0.2,
  alpha = 0.5,
  configuration = config,
  resource_landscape = fq,
  nesting_landscape = NULL,
  distance_decay = "exponential"
)

# Save results
saveRDS(result_ideal, "jenna/sim_data/idealdata.rds")
#result_ideal = readRDS("jenna/sim_data/idealdata.rds")
saveRDS(result_real, "jenna/sim_data/realdata.rds")
#result_real = readRDS("jenna/sim_data/realdata.rds")

# Write outputs to variables
yobs = result_ideal[[1]]
colony_data = result_ideal[[2]]
trap_data = result_ideal[[3]]

# # Remove undetected colonies
yobs_detected = yobs[rowSums(yobs) >0,]
colony_data_detected = colony_data[rowSums(yobs) > 0,]

# Prep data list for Stan
data = list()
data$y = yobs_detected
data$C = nrow(data$y)
data$K = ncol(data$y)
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$lowerbound = 400
data$upperbound = 1100
data$landscape = colony_data_detected$landscape_metric
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 3.5
data$rho_sd = 0.5

# Run model
stanfile = "jenna/time_and_space/multinomial_landscape.stan"

stanFit = stan(file = stanfile,
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               iter = 2000,
               warmup = 1000,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               verbose = TRUE)


########################
# Repeat with more realistic data
########################
# Write outputs to variables
yobs = result_real[[1]]
colony_data = result_real[[2]]
trap_data = result_real[[3]]

# # Remove undetected colonies
yobs_detected = yobs[rowSums(yobs) >0,]
colony_data_detected = colony_data[rowSums(yobs) > 0,]

# Prep data list for Stan
data = list()
data$y = yobs_detected
data$C = nrow(data$y)
data$K = ncol(data$y)
data$trap = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y))
data$lowerbound = 400
data$upperbound = 1100
data$landscape = colony_data_detected$landscape_metric
data$floral = trap_data$fq
data$priorVa = 1
data$priorCo = 1
data$rho_center = 3.5
data$rho_sd = 0.5

# Run model
stanfile = "jenna/time_and_space/multinomial_landscape.stan"

stanFit = stan(file = stanfile,
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               iter = 10000,
               warmup = 1000,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               verbose = TRUE)

# # for running in cmdstanr
# # add grainsize to data list
# threads_per_chain = 4
# grainsize <- max(floor(data$C / (threads_per_chain * 5)), 1)
# data$grainsize = grainsize
# 
# # compile model
# mod <- cmdstan_model(
#   "/home/melanson/projects/def-ckremen/melanson/fv_landscapeforaging/models/rs_landscape_colonylevel.stan",
#   force_recompile = TRUE, cpp_options = list(stan_threads = TRUE)
# )
# 
# #fit and save model
# print("Starting sampling.")
# fit <- mod$sample(
#   data = data,
#   chains = 4,
#   parallel_chains = 4,
#   threads_per_chain = threads_per_chain,
#   refresh = 100,
#   iter_warmup = 1000,
#   iter_sampling = 5000,
#   init = 1
# )

