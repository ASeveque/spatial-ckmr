---
title: Data, code and manuscript for 'Spatial close-kin mark-recapture models applied to terrestrial species with continuous natal dispersal'
--- 

---

* **Step-by-step guide**

  * **data:** some files required to run the mock Michigan example.
  
  * **functions:**

    * *pairwise_comparison.R:* contains two functions to create the pairwise comparison dataframes that will be used as data in the CKMR models. The first function `build.pairwise.spatial.mop` trims the set of comparisons to fit the conditional expectation of kinship. The second function `build.grouped.spatial.mop` groups mother-offspring pairs and unrelated pairs by identical covariates (essentially going from bernoulli trial of pairs to binomial distribution).
    
    * *simulate_population.R:* operate the mock Michigan population simulation and save the results as .RData.
    
    * *spatial_parameters:* contains the two crème de la crème functions to run spatial CKMR. The first function `pij` estimates natal range probability for each pairwise comparisons based on observed distances between verified mother-offspring pairs. The second function `kj` calculates the sampled area (kj) and abundance-weighted sampled area (kjD) around each potential offspring, used to derive total abundance from local abundance.
    
  * **Spatial_CKMR_guide.pdf:** guides the user through the spatial CKMR procedure. This is a great start to understand how the functions work.

* **Publication**

  * **Code**

    * *1_load_parameters.R:* load packages, source functions and population parameters required to simulate the populations and/or run CKMR models.
    * *2_simulate_populations.R:* operate a population simulation and save the results as .RData. The simulation will follow the spatial point process x sampling scenarios defined in line 19 (see source function "population_simulation.R").
    * *3_run_CKMR.R:* run the CKMR models with Nimble. You can either load an already-simulated population, or create a new one using scripts 1 and 2. It is important to run this script (#3) sequentially, as functions build on top of each other (eg the "kj" function cannot be executed if the function "pij" was not executed first).
    * *population_simulation.R:* contains six functions for the population simulations. Each scenario of spatial point process x sampling has its own function.
    * *density_maps.R:* contains two functions to create the relative abundance maps for the whole population and for the sampled individuals (to see if sampling was spatially biased).
    * *pairwise_comparison.R:* contains two functions to create the pairwise comparison dataframes that will be used as data in the CKMR models. The first function `build.pairwise.spatial.mop` trims the set of comparisons to fit the conditional expectation of kinship. The second function `build.grouped.spatial.mop` groups mother-offspring pairs and unrelated pairs by identical covariates (essentially going from bernoulli trial of pairs to binomial distribution).
    * *spatial_parameters:* contains the two crème de la crème functions to run spatial CKMR. The first function "pij" estimates natal range probability for each pairwise comparisons based on observed distances between verified mother-offspring pairs. The second function calculates the sampled area (kj) and abundance-weighted sampled area (kjD) around each potential offspring, used to derive total abundance from local abundance.
    
  * **Simulated_pops:** contains the output of the population simulations (.RData) for each of the scenarios investigated in the study.
  
---

**Abstract:** Close-kin mark–recapture (CKMR) methods use information on genetic relatedness among individuals to estimate demographic parameters. An individual’s genotype can be considered a “recapture” of each of its parent’s genotype, and the frequency of kin-pair matches detected in a population sample can directly inform estimates of abundance. CKMR inference procedures require analysts to define kinship probabilities in functional forms, which inevitably involve simplifying assumptions. Among others, population structure can have a strong influence on how kinship probabilities are formulated. Many terrestrial species are philopatric or face barriers to dispersal, and not accounting for dispersal limitation in kinship probabilities can create substantial bias if sampling is also spatially structured (e.g., via harvest). We present a spatially explicit formulation of CKMR that corrects for incomplete mixing by incorporating natal dispersal distances and spatial distribution of individuals into the kinship probabilities. We used individual-based simulations to evaluate the accuracy of abundance estimates obtained with one spatially naïve and two spatially explicit CKMR models across six scenarios with distinct spatial patterns of relative abundance and sampling probability. Estimates of abundance obtained with a CKMR model naïve to spatial structure were negatively biased when sampling was spatially biased. Incorporating patterns of natal dispersal in the kinship probabilities helped address this bias, but estimates were not always accurate depending on the model used and scenario considered. Incorporating natal dispersal into spatially structured CKMR models can address the bias created by population structure and heterogeneous sampling, but will often require additional assumptions and auxiliary data (e.g., relative abundance indices). The models shown here were designed for terrestrial species with continuous patterns of natal dispersal and high year-to-year site fidelity, but could be extended to other species.

