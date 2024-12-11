####---------- DISPERSAL AND Pij PARAMETERS ----------####

pij <- function(pop_mom_pairs) {
  
  ####----------├ Import data ----------####
  
  # Import verified mother-offspring pairs
  MomDau_true <- pop_mom_pairs[pop_mom_pairs$KinPair == "MomDau",] # Mom-Daughter
  MomSon_true <- pop_mom_pairs[pop_mom_pairs$KinPair == "MomSon",] # Mom-Son
  
  ###----------├ Create hazard rate function ----------####
  
  # Fitting hazard-rate function using nimbleDistance
  
  disp.max.F <- max(MomDau_true$distance)*3
  disp.max.M <- max(MomSon_true$distance)*3
  
  # Females
  # Model code
  y<-MomDau_true$distance
  nind<-length(y)
  distCodeV<-nimbleCode({
    y[1:nind]~dHR_V(b, sigma, Xmax = disp.max.F, point=0)
    sigma~dunif(1, disp.max.F) # sigma = standard deviation of the half-normal distribution
    b~dunif(1, 50) # b = shape of the hazard rate function
  })
  # Inits and monitors
  inits <- function() list(sigma=(disp.max.F/2), b=5)
  params <- c("sigma", "b")
  
  # Run model
  samples_F <- nimbleMCMC(
    code = distCodeV,
    constants = list(nind=nind, disp.max.F = disp.max.F),
    data = list(y=y),
    inits = inits,
    monitors = params,
    niter = 50000,
    nburnin = 10000)
 
  # Males
  # Model code
  y<-MomSon_true$distance
  nind<-length(y)
  distCodeV<-nimbleCode({
    y[1:nind]~dHR_V(b, sigma, Xmax = disp.max.M, point=0)
    sigma~dunif(1, disp.max.M)
    b~dunif(1, 50)
  })
  #inits and monitors
  inits <- function() list(sigma=(disp.max.M/2), b=5)
  params <- c("sigma", "b")
  
  samples_M <- nimbleMCMC(
    code = distCodeV,
    constants = list(nind=nind, disp.max.M = disp.max.M),
    data = list(y=y),
    inits = inits,
    monitors = params,
    niter = 50000,
    nburnin = 10000)
  
   # Extract hazard rate function parameters
  dHR_F_sigma <- mean(samples_F[,"sigma"])
  dHR_F_b <- mean(samples_F[,"b"])
  dHR_M_sigma <- mean(samples_M[,"sigma"])
  dHR_M_b <- mean(samples_M[,"b"])
  
  ###----------├ Fill pairwise comparisons dataframe  ----------####
  
  # Let's bin everything by 0.01. This reduces the number of comparisons whilst keeping pij pretty similar. 
  # This really helps binning pairwise comparisons in the CKMR model later, and making it run faster  
  
  pij_mom_daughter <- vector() # Create placeholders for pij values
  pij_mom_son <- vector()
  
  for (i in 1:nrow(mom_daughter)){ #For each Mom x Daughter comparison
    
    # The scale parameter sigma quantifies the rate of decline in kinship probability with distance between sampling locations of i and j
    pij_mom_daughter[i] <- round(dHR(mom_daughter[i,"distance"],
                                     dHR_F_b, dHR_F_sigma, Xmax = disp.max.F, point = 0, log = 0) / dHR(0, dHR_F_b, dHR_F_sigma, Xmax = disp.max.F, point = 0, log = 0),
                                 digits = 2)
    if (pij_mom_daughter[i] == 0) {pij_mom_daughter[i] <- 0.001} # Since pij is a factor in the kinship probability, it cannot be 0. 
  }
  
  for (i in 1:nrow(mom_son)){
    pij_mom_son[i] <- round(dHR(mom_son[i,"distance"],
                                dHR_M_b, dHR_M_sigma, Xmax = disp.max.M, point = 0, log = 0) / dHR(0, dHR_M_b, dHR_M_sigma, Xmax = disp.max.M, point = 0, log = 0),
                            digits = 2)
    if (pij_mom_son[i] == 0) {pij_mom_son[i] <- 0.001}
  }
  
  
  mom_daughter$pij <- pij_mom_daughter # Add the pij values to the data frame
  mom_son$pij <- pij_mom_son
  
  pop_mom_pairs <- rbind.data.frame(mom_daughter, mom_son) # Bind mother data frame
  
  return(list(pop_mom_pairs,
              dHR_F_sigma, dHR_F_b,
              dHR_M_sigma, dHR_M_b, 
              disp.max.F, disp.max.M))
}


kj <- function(density.pop.raster) {
  
  ####----------├ Prepare dataset ----------####
  
  ####----------├├ Ind 2 ----------####
  
  # Get unique offspring (drastically reduces the number of sampled_area you have to calculate). 
  # Again, M and F separated because they might have different dispersal kernels
  unique_Ind2_F <- mom_daughter[!duplicated(mom_daughter$Ind_2), c("Ind_2", "Ind_2_lon", "Ind_2_lat")]
  unique_Ind2_M <- mom_son[!duplicated(mom_son$Ind_2), c("Ind_2", "Ind_2_lon", "Ind_2_lat")]
  
  # Create a fake data set for distribution of pij with distance
  probabilities_F <- (dHR(0:disp.max.F, dHR_F_b, dHR_F_sigma, point=0) / # from 0 to max dispersal, with 1km increments
                        dHR(0, dHR_F_b, dHR_F_sigma, point=0))
  
  probabilities_M <- (dHR(0:disp.max.M, dHR_M_b, dHR_M_sigma, point=0) /
                        dHR(0, dHR_M_b, dHR_M_sigma, point=0))
  
  # Add distances to the distribution of pij
  distances_pij_F <- data.frame(distance = c(0:disp.max.F), 
                                probabilities = probabilities_F)
  distances_pij_M <- data.frame(distance = c(0:disp.max.M),
                                probabilities = probabilities_M)
  
  plot(x = distances_pij_F$distance,
       y = distances_pij_F$probabilities)
  
  # Create bins, similar to the pij bins used in pop_mom_pairs (i.e., two digits)
  # To start, you can get anything under 0.01 as a unique value, and same for 0.99
  # This will help the function run faster (cuts time in half) by removing tail of distribution
  
  # Get the distance until which pij = 0 (0.001) and 1
  # min value = max pij
  min_F <- max(distances_pij_F$distance[distances_pij_F$probabilities > 0.995]) # 0.995 because the values are rounded to 2 decimals in the hazard rate function for pij (so 0.994 is actually 0.99)
  max_F <- min(distances_pij_F$distance[distances_pij_F$probabilities < 0.005]) # 0.005 because the values are rounded to 2 decimals in the hazard rate function for pij (so 0.007 is actually 0.01)
  min_M <- max(distances_pij_M$distance [distances_pij_M$probabilities > 0.995])
  max_M <- min(distances_pij_M$distance [distances_pij_M$probabilities < 0.005]) 
  
  distances_pij_F <- distances_pij_F[min_F:max_F,]
  distances_pij_M <- distances_pij_M[min_M:max_M,]
  
  # IMPORTANT NOTE: 
  # For the simulations of the manuscript, I saved the pairwise distances as km. 
  # In the step-by-step guide, these distances are in meters. Hence why some parts of the functions are different. Using meters is probably preferable. 
  
  ####----------├├ Ind 1 ----------####
  
  # Get unique parents (just moms in this case, but real applications can look at moms and dads)
  unique_Ind1_mom <- pop_mom_pairs[!duplicated(pop_mom_pairs$Ind_1), c("Ind_1", "Ind_1_lon", "Ind_1_lat")]
  
  ####----------├ Study area ----------####
  
  # Import or create study area
  box <- matrix(c(0, 0, # Create the bounding points
                  0, 1,
                  1, 1,
                  1, 0,
                  0, 0),
                ncol=2, byrow=TRUE)
  
  sf_box <- st_polygon(list(box)) %>% # Join the points to create a polygon
    st_sfc(crs = 4326) %>% # Transform polygon into a spatial feature
    st_transform(32631) # Transform from WGS84 to UTM zone 31N (so that distances can then be calculated in km and not degrees)
  
  # Create a multipoint of your adult/mature samples (Ind 1)
  sf_unique_Ind1_mom <- st_multipoint(x = data.matrix(unique_Ind1_mom[,2:3])) %>%
    st_sfc(crs = 4326) %>%
    st_transform(32631)
  
  ####----------├ Calculate alphaj and kj ----------####
  
  ####----------├├ Daughters ----------####
  
  # Note: crs = EPSG:32631 is UTM zone 31N, for which bottom left corner is Long = 0, Lat = 0
  
  # Placeholders for scaling factor
  kj_F <- vector()
  kj_D_F <- vector()
  
  for (i in 1:nrow(unique_Ind2_F)) { # For each unique Ind2 that is female  
    
    # Create a sf point based on sampling location
    sf_Ind2_F <- st_point(x = c(unique_Ind2_F[i,"Ind_2_lon"], unique_Ind2_F[i,"Ind_2_lat"])) %>%
      st_sfc(crs = 4326) %>%
      st_transform(32631)
    
    # Placeholders for multiple spatial features for each point
    list.buffer <- list()
    list.doughnut.sf <- list()
    list.doughnut.sv <- list()
    vector.alphas <- vector()
    vector.density <- vector()
    
    # Create the buffers
    for (j in 1:nrow(distances_pij_F)) { # For each km in the fake distance data set
      
      # Create a buffer with radius = distance around the point
      list.buffer[[j]] <- st_buffer(sf_Ind2_F, dist = distances_pij_F[j,"distance"]*1000) %>% # Distance should be in meters
        st_sfc(crs = 32631)
    }
    
    # Create the doughnuts
    for (k in 2:(nrow(distances_pij_F))) { # For each of the buffers previously created
      
      # You need to manually add the first buffer (for max distance with pij = 1)
      list.doughnut.sf[[1]] <- list.buffer[[1]]  %>%
        st_intersection(sf_box) # Remove anything outside of the study area
      
      # For all other distances, remove the previous buffer (e.g., if buffer has radius 10km, remove anything under 9km radius)
      # So basically you create a doughnut with a 1km width for each distance
      list.doughnut.sf[[k]] <- st_difference(list.buffer[[k]], list.buffer[[k-1]]) %>% # Remove the section that overlaps with the inner buffers (creates a doughnut)
        st_intersection(sf_box)
      
    }
    
    # Calculate the doughnuts' areas 
    for (l in 1:length(list.doughnut.sf)) { # For each doughnut
      
      if (length(as.numeric(st_area (list.doughnut.sf[[l]]))) == 0) {break} # Will break the loop as soon as we reach a doughnut that is entirely out of the box (to avoid NA and error)
      # Note: this is only possible because max distance are hypothetical right now - not real observed dispersal values
      # Note 2: as.numeric to remove the units m^2
      
      # Get areas of each doughnut factored by pij
      # If pij = 0.5, what we are essentially saying is that "we count 50% of the adult population at this distance"
      vector.alphas[l] <- as.numeric(st_area (list.doughnut.sf[[l]])) * distances_pij_F[l, "probabilities"] 
      
    }
    
    # Calculate the mean population density (mu) for each doughnut
    for (m in 1:length(vector.alphas)) {
      
      # Convert the sf objects to SpatVector for compatiblity with terra
      list.doughnut.sv[[m]] <-  vect(list.doughnut.sf[[m]])
      
      # Get average raster cell value for each doughnut
      density <- terra::extract(density.pop.raster, list.doughnut.sv[[m]], fun = mean)
      
      # Save cell value
      vector.density[m] <- density[[2]]
      
    }
    
    # Calculate the sampled area alphaj weighted by pij only
    # Rounding helps with grouping values when doing pairwise comparisons (faster computation time)
    area_F <- sum(vector.alphas) / 1e6  # Get total "sampled" area for one point (in km^2)
    S <- st_area(sf_box) / 1e6 # Study area
    kj_F[i] <- round((as.numeric(S) / as.numeric(area_F)), digits =  2) # Ratio study area / "sampled" area
    
    # Calculate the sampled area alphaj weighted by pij AND density mu
    area_F <- sum(vector.alphas * vector.density) / 1e6  # Get total "sampled" area for one point (in km^2) after factoring each vector.area with its corresponding population density value
    S <- st_area(sf_box) / 1e6 # Study area
    kj_D_F[i] <- round((as.numeric(S) / as.numeric(area_F)), digits = 2) # Ratio study area / "sampled" area weighted by population density
    
    # Plot
    # ggplot() +
    # geom_sf(data = sf_box, fill = 'white', lwd = 0.05) +
    # geom_sf(data = list.doughnut.sf[[1]], fill = 'grey', lwd = 0.05) +
    # geom_sf(data = list.doughnut.sf[[5]], fill = 'grey', lwd = 0.05) +
    # geom_sf(data = list.doughnut.sf[[10]], fill = 'grey', lwd = 0.05) +
    # geom_sf(data = sf_Ind2_F, col = 'red', size = 1)
    
  } # Break loop over unique Ind2 (potential daughters)
  
  ####----------├├ Sons ----------####
  
  # Placeholders for scaling factor
  kj_M <- vector()
  kj_D_M <- vector()

  for (i in 1:nrow(unique_Ind2_M)) { # For each unique Ind2 that is female

    # Create a sf point based on sampling location
    sf_Ind2_M <- st_point(x = c(unique_Ind2_M[i,"Ind_2_lon"], unique_Ind2_M[i,"Ind_2_lat"])) %>%
      st_sfc(crs = 4326) %>%
      st_transform(32631)

    # Placeholders for multiple spatial features for each point
    list.buffer <- list()
    list.doughnut.sf <- list()
    list.doughnut.sv <- list()
    vector.alphas <- vector()
    vector.density <- vector()

    # Create the buffers
    for (j in 1:nrow(distances_pij_M)) { # For each km in the fake distance data set

      # Create a buffer with radius = distance around the point
      list.buffer[[j]] <- st_buffer(sf_Ind2_M, dist = distances_pij_M[j,"distance"]*1000) %>% # Distance should be in meters
        st_sfc(crs = 32631)
    }

    # Create the doughnuts
    for (k in 2:(nrow(distances_pij_M))) { # For each of the buffers previously created

      # You need to manually add the first buffer (for max distance with pij = 1)
      list.doughnut.sf[[1]] <- list.buffer[[1]]  %>%
        st_intersection(sf_box) # Remove anything outside of the study area

      # For all other distances, remove the previous buffer (e.g., if buffer has radius 10km, remove anything under 9km radius)
      # So basically you create a doughnut with a 1km width for each distance
      list.doughnut.sf[[k]] <- st_difference(list.buffer[[k]], list.buffer[[k-1]]) %>% # Remove the section that overlaps with the inner buffers (creates a doughnut)
        st_intersection(sf_box)

    }

    # Calculate the doughnuts' areas
    for (l in 1:length(list.doughnut.sf)) { # For each doughnut

      if (length(as.numeric(st_area (list.doughnut.sf[[l]]))) == 0) {break} # Will break the loop as soon as we reach a doughnut that is entirely out of the box (to avoid NA and error)
      # Note: this is only possible because max distance are hypothetical right now - not real observed dispersal values
      # Note 2: as.numeric to remove the units m^2

      # Get areas of each doughnut factored by pij
      # If pij = 0.5, what we are essentially saying is that "we counted only 50% of the adult population at this distance"
      vector.alphas[l] <- as.numeric(st_area (list.doughnut.sf[[l]])) * distances_pij_M[l, "probabilities"]

    }

    # Calculate the mean population density for each doughnut
    for (m in 1:length(vector.alphas)) {

      # Convert the sf objects to SpatVector for compatiblity with terra
      list.doughnut.sv[[m]] <-  vect(list.doughnut.sf[[m]])

      # Get average raster cell value for each doughnut
      density <- terra::extract(density.pop.raster, list.doughnut.sv[[m]], fun = mean)

      # Save cell value
      vector.density[m] <- density[[2]]

    }

    # Calculate the sampled area alphaj weighted by pij only
    area_M <- sum(vector.alphas) / 1e6  # Get total "sampled" area for one point (in km^2)
    S <- st_area(sf_box) / 1e6 # Study area
    kj_M[i] <- round((as.numeric(S) / as.numeric(area_M)), digits = 2) # Ratio study area / "sampled" area

    # Calculate the sampled area alphaj weighted by pij AND density mu
    area_M <- sum(vector.alphas * vector.density) / 1e6  # Get total "sampled" area for one point (in km^2) after factoring each vector.area with its corresponding population density value
    S <- st_area(sf_box) / 1e6 # Study area
    kj_D_M[i] <- round((as.numeric(S) / as.numeric(area_M)), digits = 2) # Ratio study area / "sampled" area weighted by population density

  } # Break loop over unique Ind2 (potential sons)
  
  ####----------├ Fill pairwise comparisons dataframe ----------####
  
  # Bring scaling factor to unique_Ind2 dataframe
  # Without relative abundance mu
  unique_Ind2_F$kj <- kj_F
  unique_Ind2_M$kj <- kj_M
  # With relative abundance mu
  unique_Ind2_F$kj_D <- kj_D_F
  unique_Ind2_M$kj_D <- kj_D_M
  
  # Join them in a unique dataframe
  unique_Ind2_join <- rbind.data.frame (unique_Ind2_F, unique_Ind2_M)

  # Reassign each scaling factor value to all pairwise comparisons (with duplicates of Ind2)
  pop_mom_pairs <- merge(pop_mom_pairs, unique_Ind2_join[,c("Ind_2", "kj", "kj_D")], by = "Ind_2")

  return(list(pop_mom_pairs))
  
}
