# Don't forget to load the population parameters and functions!

# For maximum efficiency, we run one model for each simulation with the 4 different likelihoods.
# Nf1 is the spatially-naive model; Nf2 is the natal range probability model,
# Nf3 is the natal range probability and sampled area model; Nf4 is the natal range probability and abundance-weighted sampled area model.

# We also run four chains in parallel, using a 4-core cluster

iterations <- 5 # Replication for each scenario

# Create placeholder for results
ckmr_results <- data.frame()

for(iter in 1:5) {
  set.seed(rseeds[iter])

  # I don't know why, but this function (and only this function) unloads at the end of each iter
  source("./Functions/spatial_parameters.R")
  
  ###----- Import samples and population size files -----###

  pop <- iter.loopy.list[[iter]] # List of spatial populations (only saved the population for year of estimates, to gain space)
  samples <- iter.sample.list [[iter]] # List of sampled individuals for each year of simulation
  pop.size <- iter.popsize.list [[iter]] # Population parameters for each year of simulation

  # Restrict sampling to 3 consecutive years (I saved the original files with 5 by accident)
  samples <- samples[samples$sampling.year > 27,]
  
  ###----- Create density maps -----###

  # Create a relative density map of the population.
  # This is used just for the simulations, where the entire pop is known.
  # In real life applications, one would expect to load a relative density map or some proxy
  density.pop.raster <- pop.density (pop.sf = pop)
  plot(density.pop.raster) # Look at your population
  
  # Create a relative density map of the samples
  density.samples.raster <- relative.density.samples (samples = samples)
  plot(density.samples.raster) # Look at your samples

  ###----- Create pairwise comparisons -----###

  # Run build.pairwise function, that will create HSP/POP pairwise matrices.
  pairwise <- build.pairwise.spatial (samples = samples)

  # Separate list elements into different dataframes
  pop_mom_pairs <- pairwise[[1]]

  # Separate pop comparisons by offspring sex
  mom_daughter <- pop_mom_pairs[pop_mom_pairs$Ind_2_sex == "F",]
  mom_son <- pop_mom_pairs[pop_mom_pairs$Ind_2_sex == "M",]

  # For curiosity's sake: check how many kin pairs
  MomDau <- as.numeric(count(pop_mom_pairs[pop_mom_pairs$KinPair == "MomDau",])) # Mom-Daughter
  MomSon <- as.numeric(count(pop_mom_pairs[pop_mom_pairs$KinPair == "MomSon",])) # Mom-Son

  MatPOP <- sum(MomDau, MomSon) # Mother-Offspring

  ###----- Dispersal parameters of the sample -----###
  
  ## pij ##

  # We fit a hazard-rate distance detection functions to the observed dispersal distances for verified mother-offspring pairs
  # to get the probability that sampling location of mother = natal home range offspring
  # It's an equivalent of having a scale parameter for detection probability in traditional SCR

  # Takes a minute to run
  dispersal <- pij (pop_mom_pairs = pop_mom_pairs)

  # Retrieve list elements with now pij added
  pop_mom_pairs <- dispersal[[1]]

  # Retrieve parameters for the hazard rate function for pij (needed to run next function)
  dHR_F_sigma <- dispersal[[2]]
  dHR_F_b <- dispersal[[3]]
  dHR_M_sigma <- dispersal[[4]]
  dHR_M_b <- dispersal[[5]]
  disp.max.F <- dispersal[[6]]
  disp.max.M <- dispersal[[7]]

  # Compare histograms of observed distances between all pairs (NULL) and kin pairs to detect lack of mixing
  ggplot(data = pop_mom_pairs, aes (x = distance, group = KinPair, fill = KinPair)) +
    geom_density(adjust=1.5, alpha=.4)

  ###----- Auxiliary data: spatial parameters of study site -----###

  ## Scaling factor ##
  
  # The function below calculates the sampled area and a scaling factor for each Ind2 (potential offspring) based on pij and location on the map (and relative local abundance for one of them)
  # The scaling factor is just a ratio between the state space and the sampled area. We calculate it here instead of adding it to the kinship probability, 
  # but the final result is the same (and it helps things get faster by grouping comparisons with similar scaling factor).
  
  # Note of importance: we have to be really careful with the projections used here
  # Also, this is a bit slow and probably inefficient. It could definitely be optimized.
  scaling.factor <- kj(density.pop.raster = density.pop.raster)

  # Retrieve pairwise comparisons with now kj and kj,D added
  pop_mom_pairs <- scaling.factor[[1]]
  
  ###----- Create pairwise comparison matrices -----###

  # We can't really do arrays here, since there are 5 dimensions to group. Possible, but messy.
  
  # Instead, these are pairwise comparisons (binomial distribution) grouped with identical values for the combination of
  # Ind 1 birth, Ind 2 birth, pij, kj, and kj_D
  
  # Yes is the number of MOPs; all is the total number of pairwise comparisons
 
  # We could do every single pairwise comparison as a Bernoulli trial, but the idea is to bin as much as possible for computing efficiency.
  # Here, with 5 dimensions, we still reduce the total number of comparisons by more than half.

  matrix <- build.matrix.spatial()

  pop_mom_comps <- matrix[[1]]

  ###---We write the model---###

  # Create the core cluster
  this_cluster <- makeCluster(4) # 1 MCMC = 1 core

  # Data bundle
  my.data <- list(MatPOP1 = pop_mom_comps$yes, # Positive maternal parent offspring
                  pop_mom_all_comps1 = pop_mom_comps$all, # Total maternal pop comparisons
                  MatPOP2 = pop_mom_comps$yes,
                  pop_mom_all_comps2 = pop_mom_comps$all,
                  MatPOP3 = pop_mom_comps$yes,
                  pop_mom_all_comps3 = pop_mom_comps$all,
                  MatPOP4 = pop_mom_comps$yes,
                  pop_mom_all_comps4 = pop_mom_comps$all)

  # Constants bundle
  my.constants <- list(pop_mom_length = nrow(pop_mom_comps), # Number of cohort comparisons to loop over
                       pop_mom_offspring_birthyear = pop_mom_comps$Ind_2_birth,
                       estimation_year = estimation_year,
                       pij_mom = pop_mom_comps$pij,
                       kj_mom = pop_mom_comps$kj,
                       kj_D_mom = pop_mom_comps$kj_D)

  # Model
  myCode <- nimbleCode({

    # priors
    Nf1 ~ dunif(1, 1e4) # Uniform prior for female abundance
    Nf2 ~ dunif(1, 1e4)
    Nf3 ~ dunif(1, 1e4)
    Nf4 ~ dunif(1, 1e4)

    r1 ~ dunif(-0.5, 0.5) # Uniform prior for population growth.
    r2 ~ dunif(-0.5, 0.5)
    r3 ~ dunif(-0.5, 0.5)
    r4 ~ dunif(-0.5, 0.5)

    # likelihoods
    
    # Naive model
    for (i in 1:pop_mom_length){
      MatPOP1[i] ~ dbinom(1 / (Nf1 * exp(r1 * (pop_mom_offspring_birthyear[i] - estimation_year))), pop_mom_all_comps1[i])
    }

    # Without N total derivation
    for (j in 1:pop_mom_length){
      MatPOP2[j] ~ dbinom(pij_mom[j]  / (Nf2 * exp(r2 * (pop_mom_offspring_birthyear[j] - estimation_year))), pop_mom_all_comps2[j])
    }

    # N total derivation without relative density
    for (k in 1:pop_mom_length){
      MatPOP3[k] ~ dbinom(pij_mom[k] / (Nfxj3[k] * exp(r3 * (pop_mom_offspring_birthyear[k] - estimation_year))), pop_mom_all_comps3[k])
      Nfxj3[k] <- Nf3 / kj_mom[k]
    }

    # N total derivation with relative density
    for (l in 1:pop_mom_length){
      MatPOP4[l] ~ dbinom(pij_mom[l] / (Nfxj4[l] * exp(r4 * (pop_mom_offspring_birthyear[l] - estimation_year))), pop_mom_all_comps4[l])
      Nfxj4[l] <- Nf4 / kj_D_mom[l]
    }
    
  })


  # Create a function with all the needed code
  run_MCMC_allcode <- function(seed, data, constants, code) {

    library(nimble) # Not sure why, but this needs to be specified in the function

    myModel <- nimbleModel(code = code,
                           data = data,
                           constants = constants,
                           inits = list(Nf1 = runif(1, 100, 1e4),
                                        Nf2 = runif(1, 100, 1e4),
                                        Nf3 = runif(1, 100, 1e4),
                                        Nf4 = runif(1, 100, 1e4),
                                        r1 = rnorm(1, mean = 0, sd = 0.05),
                                        r2 = rnorm(1, mean = 0, sd = 0.05),
                                        r3 = rnorm(1, mean = 0, sd = 0.05),
                                        r4 = rnorm(1, mean = 0, sd = 0.05))) # Don't need a function for inits since I will be doing 4 draws via parallel

    n.iter <- 50000
    n.burnin <- 10000
    n.thin <- 10

    CmyModel <- compileNimble(myModel)
    myMCMC <- buildMCMC(CmyModel)
    CmyMCMC <- compileNimble(myMCMC)

    results <- runMCMC(CmyMCMC, niter = n.iter, nburnin = n.burnin, thin = n.thin, setSeed = seed)

    return(results)
  }

  # Run the chains
  chain_output <- parLapply(cl = this_cluster, X = 1:4, # x = 1:n.chains
                            fun = run_MCMC_allcode,
                            data = my.data,
                            constants = my.constants,
                            code = myCode)

  # It's good practice to close the cluster when you're done with it.
  stopCluster(this_cluster)

  MCMCsummary <- MCMCsummary(chain_output, round = 3)

  MCMCsummary_Nf1 <- MCMCsummary[1,]
  MCMCsummary_Nf2 <- MCMCsummary[2,]
  MCMCsummary_Nf3 <- MCMCsummary[3,]
  MCMCsummary_Nf4 <- MCMCsummary[4,]

  Nf1_output <- c(chain_output[[1]][,"Nf1"], chain_output[[2]][,"Nf1"], chain_output[[3]][,"Nf1"], chain_output[[4]][,"Nf1"])
  Nf2_output <- c(chain_output[[1]][,"Nf2"], chain_output[[2]][,"Nf2"], chain_output[[3]][,"Nf2"], chain_output[[4]][,"Nf2"])
  Nf3_output <- c(chain_output[[1]][,"Nf3"], chain_output[[2]][,"Nf3"], chain_output[[3]][,"Nf3"], chain_output[[4]][,"Nf3"])
  Nf4_output <- c(chain_output[[1]][,"Nf4"], chain_output[[2]][,"Nf4"], chain_output[[3]][,"Nf4"], chain_output[[4]][,"Nf4"])

  iter_results <- cbind.data.frame(rseeds[iter], MatPOP,

                                   init_pop_size = init.pop.size,
                                   final_pop_size = mean(pop.size[num.years:(num.years-4),2]), # Average of the last 5 years of the simulation (pre-harvest)
                                   sample_size = nrow(samples), # Total sample size

                                   Nf_real = pop.size[estimation_year, 5],

                                   Nf1_CKMR_mode = round(density(Nf1_output)$x[which.max(density(Nf1_output)$y)], digits = 0), # Get the peak of the posterior distribution (highest density point)
                                   Nf1_CKMR_mean = MCMCsummary_Nf1$mean, Nf1_sd = MCMCsummary_Nf1$sd,
                                   `Nf1_2.5%` = MCMCsummary_Nf1$`2.5%`,`Nf1_50%` = MCMCsummary_Nf1$`50%`, `Nf1_97.5%` = MCMCsummary_Nf1$`97.5%`,
                                   Nf1_Rhat = MCMCsummary_Nf1$Rhat, Nf1_n.eff = MCMCsummary_Nf1$n.eff,

                                   Nf2_CKMR_mode = round(density(Nf2_output)$x[which.max(density(Nf2_output)$y)], digits = 0),
                                   Nf2_CKMR_mean = MCMCsummary_Nf2$mean, Nf2_sd = MCMCsummary_Nf2$sd,
                                   `Nf2_2.5%` = MCMCsummary_Nf2$`2.5%`,`Nf2_50%` = MCMCsummary_Nf2$`50%`, `Nf2_97.5%` = MCMCsummary_Nf2$`97.5%`,
                                   Nf2_Rhat = MCMCsummary_Nf2$Rhat, Nf2_n.eff = MCMCsummary_Nf2$n.eff,

                                   Nf3_CKMR_mode = round(density(Nf3_output)$x[which.max(density(Nf3_output)$y)], digits = 0),
                                   Nf3_CKMR_mean = MCMCsummary_Nf3$mean, Nf3_sd = MCMCsummary_Nf3$sd,
                                   `Nf3_2.5%` = MCMCsummary_Nf3$`2.5%`,`Nf3_50%` = MCMCsummary_Nf3$`50%`, `Nf3_97.5%` = MCMCsummary_Nf3$`97.5%`,
                                   Nf3_Rhat = MCMCsummary_Nf3$Rhat, Nf3_n.eff = MCMCsummary_Nf3$n.eff,

                                   Nf4_CKMR_mode = round(density(Nf4_output)$x[which.max(density(Nf4_output)$y)], digits = 0),
                                   Nf4_CKMR_mean = MCMCsummary_Nf4$mean, Nf4_sd = MCMCsummary_Nf4$sd,
                                   `Nf4_2.5%` = MCMCsummary_Nf4$`2.5%`,`Nf4_50%` = MCMCsummary_Nf4$`50%`, `Nf4_97.5%` = MCMCsummary_Nf4$`97.5%`,
                                   Nf4_Rhat = MCMCsummary_Nf4$Rhat, Nf4_n.eff = MCMCsummary_Nf4$n.eff)

  ckmr_results <- rbind (ckmr_results, iter_results)

  gc() # Releases memory from unused objects between loops, or R will come to an error

  MCMCtrace(chain_output, pdf = FALSE, iter = 5000) # Look at traceplots

  print(paste("Iteration", iter, "/", iterations, "finished"), sep=" ")

}  # End loop over iterations

# Calculate coefficient of variation for each iteration
ckmr_results$Nf1_CV <- (ckmr_results$Nf1_sd / ckmr_results$Nf1_CKMR_mean) * 100
ckmr_results$Nf2_CV <- (ckmr_results$Nf2_sd / ckmr_results$Nf2_CKMR_mean) * 100
ckmr_results$Nf3_CV <- (ckmr_results$Nf3_sd / ckmr_results$Nf3_CKMR_mean) * 100
ckmr_results$Nf4_CV <- (ckmr_results$Nf4_sd / ckmr_results$Nf4_CKMR_mean) * 100

# Calculate relative bias (from the mean of the distribution) for each iteration
ckmr_results$Nf1_bias <- (ckmr_results$Nf1_CKMR_mode - ckmr_results$Nf_real) / ckmr_results$Nf_real
ckmr_results$Nf2_bias <- (ckmr_results$Nf2_CKMR_mode - ckmr_results$Nf_real) / ckmr_results$Nf_real
ckmr_results$Nf3_bias <- (ckmr_results$Nf3_CKMR_mode - ckmr_results$Nf_real) / ckmr_results$Nf_real
ckmr_results$Nf4_bias <- (ckmr_results$Nf4_CKMR_mode - ckmr_results$Nf_real) / ckmr_results$Nf_real

# Save results
write.csv(ckmr_results, "ipp_restricted_2.csv", row.names = TRUE)
