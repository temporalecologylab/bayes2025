##### Incorporating multiple landscapes and timepoints into simulation and modeling workflow #####
# Script Initiated: September 15, 2025
# By: Jenna Melanson

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
library(gridExtra)
library(terra)
library(filelock)
library(tibble)

##### Set Environment #####
setwd("/Users/jenna1/Documents/UBC/bayes2025") # local
options(mc.cores = parallel::detectCores())

##### Source functions #####
source("jenna/src/GeneralizedSimFunctions.R")

##### Simulate multi-landscape, multi-timepoint data #####
result = draw_multi_landscape_timepoint(sample_size = 1000,
                                        num_landscape = 3,
                                        num_timepoints = 2,
                                        landscape_size = 1500,
                                        trapgrid_size = 300,
                                        number_traps = 25,
                                        number_colonies = 10000,
                                        colony_sizes = rep(100,10000),
                                        rho = 50,
                                        theta = 0.5,
                                        distance_decay = "exponential")
yobs_sum = result[[1]][[1]] + result[[1]][[2]]
colony_data = result[[2]]
trap_data = result[[3]]
floral = result[[4]]

saveRDS(result, "jenna/time_and_space/sim_data/observations.RDS")

# check distribution of sibship sizes
tbl <- table(rowSums(yobs_sum)[rowSums(yobs_sum) > 0])
prop_tbl <- prop.table(tbl)  # gives proportions
barplot(tbl, ylab = "Frequency (number of colonies)", xlab = "Number of siblings")

##### Fit with basic multinomial model #####
yobs = yobs_sum[rowSums(yobs_sum) >0,]

data = list()
data$y = yobs ## a C x K matrix of observations (colony C at trap K)
data$C = nrow(data$y) ## number of colonies
data$K = ncol(data$y) ## number of traps
data$lowerbound = 0 ## lower bound on spatial locations of colonies
data$upper_y = 1500 ## upper bound (in y direction)
data$upper_x = 4500 ## upper bound (in x direction)
data$floral = rowMeans(floral) ## floral resource quality
data$trap_pos = as.matrix(cbind(trap_data$trap_x, trap_data$trap_y)) ## trap coordinates


#select stan model to fit
stanfile = "jenna/time_and_space/multinomial.stan"

#fit and save model
stanFit = stan(file = stanfile,
               data = data, seed = 5838299,
               chains = 4, cores = 4,
               iter = 6000, warmup = 1000,
               verbose = TRUE)
saveRDS(stanFit, "jenna/time_and_space/sim_data/simplemodel.RDS")