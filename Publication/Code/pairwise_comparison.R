####---------- Pairwise comparisons ----------####

build.pairwise.spatial <- function(samples){
  
  # Create a matrix of distances between all samples
  # For now, linear distances. But might need to find a better way for real uses (e.g., based on connectivity).

  samples.sf <- st_as_sf (samples, coords = c("lon", "lat"), crs = 4326) # Transform your dataframe as a sf

  samples.distances.matrix <- st_distance(x = samples.sf, y = samples.sf) # Calculate distances between two individuals. Surprisingly fast with linear distances

  colnames (samples.distances.matrix) <- samples$ID # Transform col and row names to the sample names
  rownames (samples.distances.matrix) <- samples$ID

  samples.distances.matrix[1:10, 1:10] # Looks like this (too big to show the entire matrix)

  # Create dataframe of pairwise comparisons with just individual IDs

  pairwise.df <- data.frame(t(combn(samples$ID, m=2))) # Generates all the possible combinations. m = 2 means pair comparisons. t(x) gives the dimensions of the matrix (dataframe).
  colnames(pairwise.df) <- c("Ind_1", "Ind_2") # Rename columns so they can easily be joined
  head(pairwise.df)
  
  distance <- vector() # Create placeholder for pairwise distances
  
  for(i in 1:nrow(samples)-1){
  dist <- as.vector(samples.distances.matrix[i,(i+1):nrow(samples)]) # Extract relevant pairwise distances for each Ind_1. The as.vector removes the unit automatically.
  distance <- append(distance, dist)
  }
  
  # Assign distances to the pairwise df
  pairwise.df$distance <- distance

  # For diagnostics: compare the distance matrix and the dataframe (the two last values should match)
  test <- pairwise.df[sample(1:nrow(pairwise.df), 1),] # Select a random row of pairwise df
  test$distance
  samples.distances.matrix[rownames(samples.distances.matrix) == test[1,1], colnames(samples.distances.matrix) == test[1,2]] # Independently look at distance between these two ind in distance matrix
  
  # Bring all other relevant columns for Ind 1 and Ind 2

  # Create dataframe that will be used to extract the birth years for the younger individual from each pairwise comparison using joins.
  Ind1_birthyears <- samples[,c("ID", "birth.year", "sampling.year", "mother", "father", "sex", "lon", "lat")] # Select relevant columns only
  colnames(Ind1_birthyears) <- c("Ind_1","Ind_1_birth","Ind_1_sampling", "Ind_1_mom", "Ind_1_dad", "Ind_1_sex", "Ind_1_lon", "Ind_1_lat") # Rename columns (make sure they are in the same order as the previous command)
  
  Ind2_birthyears <- samples[,c("ID", "birth.year", "sampling.year", "mother", "father", "sex", "lon", "lat")] # Select relevant columns only
  colnames(Ind2_birthyears) <- c("Ind_2","Ind_2_birth","Ind_2_sampling", "Ind_2_mom", "Ind_2_dad", "Ind_2_sex", "Ind_2_lon", "Ind_2_lat") # Rename columns (make sure they are in the same order as the previous command)
  
  # Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
  # This is the main pairwise comparison matrix with all comparisons (not all relevant for now) and individual data.
  
  pairwise.df_all <- merge(x = pairwise.df, y = Ind1_birthyears, 
                           by = "Ind_1", sort = FALSE)

  pairwise.df_all <- merge(x = pairwise.df_all, y = Ind2_birthyears, 
                           by = "Ind_2", sort = FALSE)
  
  # Change order of columns (for ease of read and logic)
  pairwise.df_all <- pairwise.df_all[, c("Ind_1", "Ind_1_sex", "Ind_1_birth", "Ind_1_sampling", "Ind_1_mom", "Ind_1_dad", "Ind_1_lon", "Ind_1_lat",
                                         "Ind_2", "Ind_2_sex", "Ind_2_birth", "Ind_2_sampling", "Ind_2_mom", "Ind_2_dad", "Ind_2_lon", "Ind_2_lat",
                                         "distance")]
  
  # Switch distances to km (easier to read)
  pairwise.df_all$distance <- pairwise.df_all$distance/1000
  
  head(pairwise.df_all)
  
  # Remove comparisons that are "not biologically possible"
  
  # Important: only adults should be sampled as potential parents (Ind 1). 
  # If juveniles are sampled (genotyped), N goes up quasi exponentially (see Waples & Feutry, 2021)
  
  pop.pairwise.df_filt <- pairwise.df_all %>%
    filter ((Ind_1_sampling - Ind_1_birth) >= repro.age) %>% # Remove potential parents that were harvested before reaching age at maturity
    filter ((Ind_2_birth - Ind_1_birth) >= repro.age) %>% # Remove potential parents that did not reach age at maturity when offspring was born
    filter ((Ind_2_birth - Ind_1_birth) <= max.age) %>% # Remove potential parents that were dead when offspring was born (can happen if more than one sampling occasion)
    filter (Ind_1_sampling >= Ind_2_birth) %>% # Remove potential parents that were sampled before offspring was born (because of lethal sampling). For now, let's assume that adults can reproduce and then die in the same year.
    filter(Ind_1_birth != Ind_2_birth) #Filter intra-cohort comparisons (redundant with line 2 but keep it for clarity's sake)
  
  # Separate mothers
  pop.pairwise.df_filt_female <- pop.pairwise.df_filt[pop.pairwise.df_filt$Ind_1_sex == "F",]

  # Assign kinships
  pop.pairwise.df_filt_female$KinPair <- ifelse(pop.pairwise.df_filt_female$Ind_1 == pop.pairwise.df_filt_female$Ind_2_mom & pop.pairwise.df_filt_female$Ind_2_sex == "F", "MomDau", # Mother-daughter
                                        ifelse(pop.pairwise.df_filt_female$Ind_1 == pop.pairwise.df_filt_female$Ind_2_mom & pop.pairwise.df_filt_female$Ind_2_sex == "M", "MomSon", # Mother-son
                                               "UR")) # Unrelated

   return(list(pop.pairwise.df_filt_female)
          
  )
}


####---------- Matrices ----------####

build.matrix.spatial <- function() {
  
  # Find positive mother-offspring pairs
  pop_mom_positives <- pop_mom_pairs %>%
    filter(Ind_1 == Ind_2_mom) %>% 
    dplyr::select(Ind_1_birth, Ind_2_birth, pij, kj, kj_D) %>% 
    plyr::count() 
  
  # Find negatives
  pop_mom_negatives <- pop_mom_pairs %>%
    filter(Ind_1 != Ind_2_mom) %>%
    dplyr::select(Ind_1_birth, Ind_2_birth, pij, kj, kj_D) %>% 
    plyr::count()
  
  # Join and group by similar comparisons
  pop_mom_comps <- pop_mom_positives %>% 
    dplyr::rename(yes = freq) %>% 
    full_join(pop_mom_negatives, by = c("Ind_1_birth", "Ind_2_birth", "pij", "kj", "kj_D")) %>% 
    dplyr::rename(no = freq) %>% 
    mutate(yes = replace_na(yes, 0), no = replace_na(no, 0)) %>%
    mutate(all = yes + no)
  
  return(list(pop_mom_comps))
  
}
