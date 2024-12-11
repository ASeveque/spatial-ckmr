####---------- Do science good ----------####

# Don't forget to load the population parameters and functions! 

iterations <- 100 # Replication for each scenario
rseeds <- sample(0001:500000, iterations) # One different seed for each iteration

####---------- Population simulation ----------####

# Create placeholders for all samples of iterations
iter.loopy.list <- list() 
iter.sample.list <- list() 
iter.popsize.list <- list()

for(iter in 1:iterations) {
  set.seed(rseeds[iter])
  
  # Run population simulation
  popsim <- HPP.gradient()
  
  # Extract relevant output dataframes
  pop <- popsim [[1]] # List of data frames for each year of simulation
  sample.list <- popsim [[2]] # List of sampled individuals for each year of simulation
  pop.size <- popsim [[3]] # Population parameters for each year of simulation
  
  # Verify that your population is stable (plus minus stochasticity)
  par(mfrow=c(2,2))
  plot(pop.size$population_size_preharvest , pch=19, ylim = c(0.9*min(pop.size$population_size_preharvest ), 1.1*max(pop.size$population_size_preharvest))) #Plot population size through time
  plot(pop.size$Nm_preharvest, pch=19, ylim = c(0.9*min(pop.size$Nm_preharvest), 1.1*max(pop.size$Nm_preharvest))) #Plot Nm through time
  plot(pop.size$Nf_preharvest, pch=19, ylim = c(0.9*min(pop.size$Nf_preharvest), 1.1*max(pop.size$Nf_preharvest))) #Plot Nm through time

  ###----- Samples selection -----###
  
  # Select your samples that will be used to run the CKMR.
  
  samples <- data.frame() # Placeholder
  
  for(s in sample_subset){ # For each year of sampling
    samples_year <- sample.list[[s]] # Retrieve sampled individuals
    samples_year$sampling_year <- rep(s) # Keep track of what year they were sampled
    samples <- rbind(samples, samples_year) # Bind
  }
  
  # Sort dataframe by birth year, so that individual 1 is always older than individual 2 (because it comes first in the df)
  samples <- samples[order(samples$birth.year),]

  ###----- Store simulation output -----###

  iter.loopy.list[[iter]] <- pop[[estimation_year]]  # Store the entire population at year of estimation (if you keep years, the RData gets quite big)
  iter.sample.list[[iter]] <- samples # Store sampled individuals
  iter.popsize.list[[iter]] <- pop.size # Store yearly population size

  gc() # Releases memory from unused objects between loops
  
  print(paste("Iteration", iter, "/", iterations, "finished"), sep=" ")
  
}

# Save simulation output as RData
save(rseeds, iter.loopy.list, iter.sample.list, iter.popsize.list, file = "IPP_uniform.RData")
