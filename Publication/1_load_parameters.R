####---------- Load packages and source functions ----------####

library(tidyverse) # Data management
library(nimble) # Run models
library(MCMCvis) # Visualize models
library(parallel) # For parallel processing of MCMCs
library(DescTools) # To get mode from posterior distribution
library(sf) # Spatial stuff
library(terra) # Spatial stuff (for rasters)
library(geosphere) # Spatial stuff (for simulations only)
library(fitdistrplus) # To obtain dispersal pdf
library(Distance) # To obtain dispersal pdf
library(nimbleDistance) # To obtain dispersal pdf
library(ggplot2) # For plots
library(gridExtra) # For plots
library(cowplot) # For plots
library(scales) # For plots
library(ggforce) # For plots

source("./Functions/population_simulation.R")
source("./Functions/pairwise_comparison.R")
source("./Functions/spatial_parameters.R")
source("./Functions/density_maps.R")

####---------- Demographic parameters ----------####

# Initial population
init.prop.female <- 0.5 # Proportion of the initial population size that is female
repro.age <- 2 # Set age of reproductive maturity
max.age <- 10 # Set the maximum age allowed in the simulation

num.years <- max.age * 3 # The number of years to run in the simulation

# Year 0 breeding
birth.sex.ratio <- c(0.5,0.5) # The probability that each baby is F:M
num.mates <- c(1:2) # Possible multiple paternity

# All breeding years
yoy.survival <- 0.5 # Mean young of year survival ish (just for Leslie matrix)
juvenile.survival <- 0.57 # Juvenile survival - between 1 years old and repro.age
adult.survival <- 0.57 # Adult survival - after age of repro.age
init.adult.pop.size <- 863 # Initial mothers population

# Create an age distribution based on average fecundity and survival
f <- 3 # Average annual adult fecundity (= number of offspring)
props <- rep (NA, max.age+1) # Placeholder
props[1] <- f # YOY (those just being born)
props[2] <- f * yoy.survival # YOY survival (from 0 to 1 years)
for (y in 3:(repro.age+1)) props[y] <- props[y-1] * juvenile.survival # Juvenile survival (from 1 to adult years)
for (y in (repro.age+2):(max.age+1)) props[y] <- props[y-1] * adult.survival # Adult survival

Nages <- round(props[-1] * init.adult.pop.size) # Transforms age "ratios" to number of individuals in each age + remove YOY
init.pop.size <- sum(Nages) # Total population size = 2,999

# Sampling variables
harvest.size = 10 # Harvest n% of the previous year's population
harvest.age = 1 # Minimum age of harvest
sample_first_y = num.years - 2 # First year from simulation where we start genetic sampling
sample_subset <- seq(from=sample_first_y, to=num.years, by=1) # Years of simulation where we will keep samples for CKMR. "by=n" is the sampling frequency.
estimation_year = 28