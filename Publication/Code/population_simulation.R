HPP.uniform <- function(){
  
  # This is an attempt to "force" population to be homogeneous (i.e., force some dispersal on the outside of the box)
  
  ###---Set up your initial population---###
  
  init.pop <- data.frame() # Create a blank data frame which will become the initial population
  
  age <- NULL # Create an empty vector for age
  for(y in 1:length(Nages)) age <- c(age, rep(y, Nages[y])) # Elongate age into population size and fill each row according to ratios in Nages
  
  # Loop that creates the below data for each individual in a population the size of "init.pop.size"
  for(i in 1:init.pop.size){
    
    indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
    mother <- "xxxxx" # The individuals in the initial population do not have known mothers
    father <- "xxxxx" # The individuals in the initial population do not have known fathers
    birth.year <- -1 # this is a place holder for individuals born within the simulation
    sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the proportions set in the parameter section
    dispersal <- 1 # Let's assume they all dispersed, no matter the age (only 1 dispersal event in life)
    init.vec <- cbind.data.frame(indv.name, birth.year, age[i], mother, father, sex, dispersal) # Create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
    init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
  }
  
  names(init.pop)[3] <- "age" # Change the name to remove the [i] from the loop
  
  # Have a look at your population
  
  head(init.pop)
  summary(as.factor(init.pop$age)) # Should match Nages
  summary(as.factor(init.pop$sex)) # Should be 50/50
  
  # Create a box (100x100km)
  box <- matrix(c(0, 0, # Create the bounding points
                  0, 1,
                  1, 1,
                  1, 0,
                  0, 0),
                ncol=2, byrow=TRUE)
  
  sf_box <- st_polygon(list(box)) # Join the points to create a polygon
  sf_box <- st_sfc(sf_box, crs = 4326) # Transform polygon into a spatial feature
  
  # Then create two grids on box
  # One larger grid for reproduction
  sf_grid_large  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.1, 0.1)) %>% # Make it 10x10
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326) # Reverts the dataframe to a sf object (now including ID)
  
  # One smaller grid for dispersal and heterogenous point process (/density)
  sf_grid_small  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.05, 0.05)) %>% # Make it 50x50
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326)
  
  # Get centroid for smaller grid
  sf_grid_small_centroid <- st_centroid(sf_grid_small)
  
  # Create random points that will populate the small grid (centroid of each cell)
  points_init <- st_sample(sf_grid_small_centroid, size = nrow(init.pop), replace = TRUE)
  
  # Assign those points to the init.pop dataframe
  init.pop$long <- st_coordinates(points_init)[,1]
  init.pop$lat <- st_coordinates(points_init)[,2]
  
  # Assign each ind to the cell they are into
  init.pop.sf <- st_as_sf (init.pop, coords = c("long", "lat"), crs = 4326) %>% # Transform your dataframe as a sf
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  # Rename columns for convenience
  names(init.pop.sf)[names(init.pop.sf) == 'ID.x'] <- 'ID.small'
  names(init.pop.sf)[names(init.pop.sf) == 'ID.y'] <- 'ID.large'
  
  # Visualize what you just did
  
  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_large, lwd = 0.3) +
    geom_sf(data = sf_grid_small, lwd = 0.3) +
    geom_sf(data = init.pop.sf, mapping = aes (color = sex), size = 0.7)
  
  # Create matrix and fill with distances
  distance.matrix <- st_distance(x = sf_grid_small_centroid, y = sf_grid_small_centroid) # units = m
  distance.matrix[1:10, 1:10] # ~5.5km between adjacent grids
  
  # Transition matrix based on distances only (follows dgamma with shape = 2 and rate = 1.5e-4)
  # Note: might need to do one matrix for males and one matrix for females if dispersal distances are different
  transition.distance = matrix(dgamma(as.numeric(distance.matrix), 2, 1.5e-4), 400, 400)
  transition.distance[1:10, 1:10]
  
  # Divide by a normalizing constant so that sum(transition.distance) == 1
  transition.distance <- transition.distance / sum(transition.distance)

  # Calculate SUM probability of ind moving to each cell based on distance only
  transition.distance.vector <- colSums(transition.distance) 
  sf_grid_small_centroid$TDV <- transition.distance.vector # Using the plot below, you can see that there will naturally be more points in the center of the box
  
  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_small_centroid, mapping = aes (color = TDV))
  
  # Create a density that is the *exact* opposite of the natural density created by distances
  # We are basically forcing a perfect density homogeneity (minus stochastic)
  mu_opposite <-  mean(transition.distance.vector) / transition.distance.vector
  
  # Transition matrix with distances factored with mu 
  transition.distance.mu <- matrix(0, 400, 400) 
  
  for (i in 1:nrow(sf_grid_small)) { # For each row ("from")
    for (j in 1:nrow(sf_grid_small)) { # For each column ("to")
      transition.distance.mu[i,j] <- transition.distance[i,j] * mu_opposite[j] # Multiply distance with mu
    }
  }
  
  transition.distance.mu[1:20, 1:20]
  
  colSums(transition.distance.mu) # Now you can see that each cell has equal probability of being "moved to".
  
  ###---Breeding---###
  
  ###---Year 0 breeding---###
  
  mothers <- which(init.pop.sf$sex=='F' & init.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
  
  YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
  
  # Loop through all the available mothers
  for(j in 1:length(mothers)){
    
    num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    
    # Record the mother
    mother <- as_Spatial(init.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
    
    # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
    fathers <- which(init.pop.sf$sex=='M' & init.pop.sf$age>=repro.age & init.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
    
    if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
    
    # Need to make a combination of number of available fathers and number of mates this year for this female
    ifelse(length(fathers) == 1,  # If only one father available...
           {num.mates.x <- 1; fathers.ID <- init.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
           ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
    
    num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
    num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
    
    if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
    
    # If there are two fathers, need to distribute the cubs among them
    ifelse (num.offspring %% 2 == 0, # If even number of cubs....
            num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
            num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
    
    for (h in 1:num.mates.x){
      
      ifelse (num.mates.x == 1, 
              num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
              ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                     num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                     num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
      
      if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
      
      indv.name <- NA # Create a place holder for the random name given to each offspring
      age <- 0 # Assign a 0 age to each offspring born this year
      birth.year <- 0 # Note that these pups were born in the 0 year of the simulation
      sex <- NA # Create a place holder for the sexes
      long <- NA
      lat <- NA
      dispersal <- NA
      
      
      inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                   "mother" = mother[[1]], # record the mother's name
                                   "father" = fathers.ID[[1]][[h]], # record the father's name,
                                   sex, long, lat, dispersal) 
      
      inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
      inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
      inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
      inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
      inter.df$dispersal <- 0 # Newborns have not dispersed yet
      
      baby.names <- vector() # Create a blank vector for random baby names
      
      for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
        name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
        baby.names <- c(baby.names,name1)
      } # End loop over name
      
      inter.df$indv.name <- baby.names # add the baby names to the data frame
      
      YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
      
    } # End loop over mates
    
  } # End loop over mothers
  
  YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
  
  # Assign each YOY to the cell they are into
  YOY.df.sf <- YOY.df.sf %>% 
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  
  # Rename columns for convenience
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'  
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
  
  
  loopy.pop.sf <- NULL
  loopy.pop.sf <- rbind(init.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
  
  # At end of year 0
  Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  
  # MORTALITY (two steps)
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
  loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                  ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                         ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
  
  # Density-dependent yoy survival
  next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
  diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
  ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
  yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
  loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
  
  # Lethal sampling of individuals. 
  
  # Size of this year's adult & juvenile population (before survival draw).
  number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
  
  sample.pop <- loopy.pop.sf %>%
    subset (survival == "M") %>% # Only harvest individuals
    subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
    slice_sample(n = round((harvest.size*number.adults)/100), replace = FALSE) # Harvest n% of this year's total adult population
  
  # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
  loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                   "H", # Then write him as harvested
                                   loopy.pop.sf$survival) # Otherwise keep S or M as it was before
  
  # Dispersal and grid assignment 
  # Since all my initial pop has dispersed, and none of the YOY for year 0 disperse, there is no dispersal event for the first year
  
  # At end of year 0
  print(paste("year 0 ", "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
  
  # Note: not all males or females are fathers or mothers
  # Mothers and fathers may have contributed to reproduction this year and die afterwards (since survival is estimated at a later stage)
  
  # Visualize the YOY
  # Don't worry if it looks like there are more males than females. The sibling dots overlap, and males are "on top" of females, but sex ratio is okay.
  # ggplot() +
  #   geom_sf(data = sf_box, fill = 'white', lwd = 0.05) +
  #   geom_sf(data = sf_grid, fill = 'transparent', lwd = 0.3) +
  #   geom_sf(data = YOY.df.sf, mapping = aes (color = sex), size = 0.7) #+
  #coord_sf(datum = st_crs(32616)) # To change axis to UTM (meters)
  
  ###---Loop for all other breeding years---###
  
  pop.size <- data.frame()
  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year
  sample.list <- list() # Make list to store dataframe of sampled individual for each year
  
  for(v in 1:(num.years)){ #loop through all of the years in the simulation
    
    v.pop.sf <- loopy.pop.sf[loopy.pop.sf$survival =="S",] # Bring in the data from the previous iteration, but only include those that survive
    v.pop.sf <- subset (v.pop.sf, select = - survival) # Remove the survival column
    v.pop.sf$age <- v.pop.sf$age+1 # Increase each individuals age by one for the new year
    
    mothers <- which(v.pop.sf$sex=='F' & v.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
    
    YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
    
    # Loop through all the available mothers
    for(j in 1:length(mothers)){
      
      num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
      inter.df <- data.frame() # Create a data frame for the offspring born from each mother
      
      # Record the mother
      mother <- as_Spatial(v.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
      
      # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
      fathers <- which(v.pop.sf$sex=='M' & v.pop.sf$age>=repro.age & v.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
      
      if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
      
      # Need to make a combination of number of available fathers and number of mates this year for this female
      ifelse(length(fathers) == 1,  # If only one father available...
             {num.mates.x <- 1; fathers.ID <- v.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
             ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
      
      num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
      num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
      
      if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
      
      # If there are two fathers, need to distribute the cubs among them
      ifelse (num.offspring %% 2 == 0, # If even number of cubs....
              num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
              num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
      
      for (h in 1:num.mates.x){
        
        ifelse (num.mates.x == 1, 
                num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
                ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                       num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                       num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
        
        if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
        
        indv.name <- NA # Create a place holder for the random name given to each offspring
        age <- 0 # Assign a 0 age to each offspring born this year
        birth.year <- v # Note that these pups were born in the 0 year of the simulation
        sex <- NA # Create a place holder for the sexes
        long <- NA
        lat <- NA
        dispersal <- NA
        
        
        inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                     "mother" = mother[[1]], # record the mother's name
                                     "father" = fathers.ID[[1]][[h]], # record the father's name,
                                     sex, long, lat, dispersal) 
        
        inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
        inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
        inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
        inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
        inter.df$dispersal <- 0 # Newborns have not dispersed yet
        
        baby.names <- vector() # Create a blank vector for random baby names
        
        for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
          name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
          baby.names <- c(baby.names,name1)
        } # End loop over name
        
        inter.df$indv.name <- baby.names # add the baby names to the data frame
        
        YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
        
      } # End loop over mates
      
    } # End loop over mothers
    
    YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
    
    # Assign each YOY to the cell they are into
    YOY.df.sf <- YOY.df.sf %>% 
      st_join(sf_grid_small, join = st_intersects)  %>% # Spatial join of each point to the grid ID. Adds a new column.
      st_join(sf_grid_large, join = st_intersects)
    
    
    # Rename columns for convenience
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
    
    loopy.pop.sf <- NULL
    loopy.pop.sf <- rbind(v.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
    
    # At end of year 0
    Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    
    # DISPERSAL
    
    # Let's fix dispersal at 1yo. 
    dispersers.pop <- data.frame() # Create a place holder for this year's dispersers
    
    dispersers.pop <- loopy.pop.sf %>%
      subset (dispersal == 0) %>% # Only one dispersal events, so remove those that have already dispersed (only used if eg only 50% of 1yo dispersed, and the rest at 2yo)
      subset (age == 1) %>% # All 1 year olds disperse
      st_drop_geometry() %>% # Remove geometry for now
      subset (select = - ID.large) # Remove large cell ID (it might change after dispersal)
    
    for (d in 1:nrow(dispersers.pop)) {
      
      dest.vec <- transition.distance.mu[as.numeric(dispersers.pop[d, "ID.small"]),] # Get all possible destinations from the cell of origin (row)
      
      dispersers.pop[d, "ID.small"] <- sample(1:nrow(sf_grid_small), size = 1, prob = dest.vec) # Select one final destination based on distance x mu probability
      
      # Retrieve geometry from small grid ID (based on centroids), transformed into long and lat columns
      dispersers.pop[d, "long"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "X"]
      dispersers.pop[d, "lat"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "Y"]
    }
    
    dispersers.pop.sf <- st_as_sf(dispersers.pop, coords = c("long", "lat"), crs = 4326) %>%  # Create the spatial feature for dispersers
      st_join(sf_grid_large, join = st_intersects) # Spatial join of each point to the large grid ID. Adds a new column.
    
    # Rename column for convenience
    names(dispersers.pop.sf)[names(dispersers.pop.sf) == 'ID'] <- 'ID.large'
    
    # Assign dispersal as being done
    dispersers.pop.sf$dispersal <- 1
    
    # # DIAGNOSTIC: let's connect the "start" and "end" of dispersal. Gives an idea of how far individuals move around.
    # 
    # dispersers.pop.origin <- loopy.pop.sf %>% # This one is just for diagnostic purposes
    #   subset (dispersal == 0) %>%
    #   subset (age == 1)
    # 
    # start.long <- st_coordinates(dispersers.pop.origin)[,1]
    # start.lat <- st_coordinates(dispersers.pop.origin)[,2]
    # 
    # end.long <- st_coordinates(dispersers.pop.sf)[,1]
    # end.lat <- st_coordinates(dispersers.pop.sf)[,2]
    # 
    # # plot(c(start.long, end.long), c(start.lat, end.lat), col = rep(c("red", "blue"), each = length(start.long)),
    # #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # # arrows(start.long, start.lat, end.long, end.lat)
    # 
    # # If it's too messy, you can subset the dataframes
    # 
    # start.long.sub <- start.long[1:(length(start.long)/5)]
    # start.lat.sub <- start.lat[1:(length(start.lat)/5)]
    # 
    # end.long.sub <- end.long[1:(length(end.long)/5)]
    # end.lat.sub <- end.lat[1:(length(end.lat)/5)]
    # 
    # plot(c(start.long.sub, end.long.sub), c(start.lat.sub, end.lat.sub), col = rep(c("red", "blue"), each = length(start.long.sub)),
    #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # arrows(start.long.sub, start.lat.sub, end.long.sub, end.lat.sub)

    # Remove individuals that are in dispersers.pop from the loopy.pop, and add them again with the new dataset
    loopy.pop.sf <- loopy.pop.sf[!loopy.pop.sf$indv.name %in% dispersers.pop.sf$indv.name,] # Remove
    loopy.pop.sf <- rbind.data.frame(loopy.pop.sf, dispersers.pop.sf) # Add
    
    # MORTALITY (two steps)
    
    #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
    loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                    ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                           ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
    
    # Density-dependent yoy survival
    next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
    diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
    ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
    yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
    loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
    
    
    
    # Lethal sampling of individuals. 
    
    # Size of this year's adult & juvenile population (before survival draw).
    number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
    
    sample.pop <- loopy.pop.sf %>%
      subset (survival == "M") %>% # Only harvest individuals
      subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
      slice_sample(n = round((harvest.size*number.adults)/100), replace = FALSE) # Harvest n% of this year's total adult population
    
    # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
    loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                     "H", # Then write him as harvested
                                     loopy.pop.sf$survival) # Otherwise keep S or M as it was before
    
    
    # Return samples to a data frame, with geometry transformed into long and lat columns
    sample.pop <- sample.pop %>%
      mutate(long = st_coordinates(.)[,"X"], #Re-extract long and lat coordinates in separate columns
             lat = st_coordinates(.)[,"Y"]) %>%
      st_set_geometry(NULL) # Remove geometry
    
    sample.list[[v]] <- sample.pop
    loopy.list[[v]] <- loopy.pop.sf # Save the current year's population data as a list element, where the index corresponds to the year
    
    # Print the simulation year and the population size in the R console so they can be observed
    print(paste("year", v, "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df.sf), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
    
    ### Count survivors for each year v
    pop.size.vec <- cbind.data.frame(year=v,
                                     population_size_preharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]), # Minus young of the year
                                     population_size_postharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1 & loopy.pop.sf$survival=="S",]), # Minus young of the year
                                     Nm_preharvest = Nm, # Includes males that may have reproduced this year but were subsequently killed (harvest or natural)
                                     Nf_preharvest = Nf) # Includes females that may have reproduced this year but were subsequently killed (harvest or natural)
    
    pop.size <- rbind(pop.size, pop.size.vec)
    
  } # end loop over sim years
  
  # Label the list elements with the year
  names(loopy.list) <- paste0("year.end.pop.", seq(1:num.years))
  
  return(invisible(list(loopy.list, sample.list, pop.size)))
  
}


HPP.gradient <- function(){
  
  # This is an attempt to "force" population to be homogeneous (i.e., force some dispersal on the outside of the box)
  
  ###---Set up your initial population---###
  
  init.pop <- data.frame() # Create a blank data frame which will become the initial population
  
  age <- NULL # Create an empty vector for age
  for(y in 1:length(Nages)) age <- c(age, rep(y, Nages[y])) # Elongate age into population size and fill each row according to ratios in Nages
  
  # Loop that creates the below data for each individual in a population the size of "init.pop.size"
  for(i in 1:init.pop.size){
    
    indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
    mother <- "xxxxx" # The individuals in the initial population do not have known mothers
    father <- "xxxxx" # The individuals in the initial population do not have known fathers
    birth.year <- -1 # this is a place holder for individuals born within the simulation
    sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the proportions set in the parameter section
    dispersal <- 1 # Let's assume they all dispersed, no matter the age (only 1 dispersal event in life)
    init.vec <- cbind.data.frame(indv.name, birth.year, age[i], mother, father, sex, dispersal) # Create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
    init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
  }
  
  names(init.pop)[3] <- "age" # Change the name to remove the [i] from the loop
  
  # Have a look at your population
  
  head(init.pop)
  summary(as.factor(init.pop$age)) # Should match Nages
  summary(as.factor(init.pop$sex)) # Should be 50/50
  
  # Create a box (100x100km)
  box <- matrix(c(0, 0, # Create the bounding points
                  0, 1,
                  1, 1,
                  1, 0,
                  0, 0),
                ncol=2, byrow=TRUE)
  
  sf_box <- st_polygon(list(box)) # Join the points to create a polygon
  sf_box <- st_sfc(sf_box, crs = 4326) # Transform polygon into a spatial feature
  
  # Then create two grids on box
  # One larger grid for reproduction
  sf_grid_large  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.1, 0.1)) %>% # Make it 10x10
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326) # Reverts the dataframe to a sf object (now including ID)
  
  # One smaller grid for dispersal and heterogenous point process (/density)
  sf_grid_small  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.05, 0.05)) %>% # Make it 50x50
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326)
  
  # Get centroid for smaller grid
  sf_grid_small_centroid <- st_centroid(sf_grid_small)
  
  # Create random points that will populate the small grid (centroid of each cell)
  points_init <- st_sample(sf_grid_small_centroid, size = nrow(init.pop), replace = TRUE)
  
  # Assign those points to the init.pop dataframe
  init.pop$long <- st_coordinates(points_init)[,1]
  init.pop$lat <- st_coordinates(points_init)[,2]
  
  # Assign each ind to the cell they are into
  init.pop.sf <- st_as_sf (init.pop, coords = c("long", "lat"), crs = 4326) %>% # Transform your dataframe as a sf
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  # Rename columns for convenience
  names(init.pop.sf)[names(init.pop.sf) == 'ID.x'] <- 'ID.small'
  names(init.pop.sf)[names(init.pop.sf) == 'ID.y'] <- 'ID.large'
  
  # Visualize what you just did
  
  sf_grid_small_centroid.coord <- as.data.frame(sf::st_coordinates(sf_grid_small_centroid))
  sf_grid_small_centroid.coord$ID <- sf_grid_small_centroid$ID
  
  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_large, lwd = 0.3) +
    geom_sf(data = sf_grid_small, lwd = 0.3) +
    #geom_sf(data = init.pop.sf, mapping = aes (color = sex), size = 0.7) +
    # geom_sf(data = sf_grid_small_centroid, size = 0.7) +
    geom_text(data = sf_grid_small_centroid.coord, aes(X, Y, label = ID), colour = "black")
  
  # Create matrix and fill with distances
  distance.matrix <- st_distance(x = sf_grid_small_centroid, y = sf_grid_small_centroid) # units = m
  distance.matrix[1:10, 1:10] # ~5.5km between adjacent grids
  
  # Transition matrix based on distances only (follows dgamma with shape = 2 and rate = 1.5e-4)
  # Note: might need to do one matrix for males and one matrix for females if dispersal distances are different
  transition.distance = matrix(dgamma(as.numeric(distance.matrix), 2, 1.5e-4), 400, 400)
  transition.distance[1:10, 1:10]
  
  # Divide by a normalizing constant so that sum(transition.distance) == 1
  transition.distance <- transition.distance / sum(transition.distance)
  
  # Calculate SUM probability of ind moving to each cell based on distance only
  transition.distance.vector <- colSums(transition.distance) 
  sf_grid_small_centroid$TDV <- transition.distance.vector # Using the plot below, you can see that there will naturally be more points in the center of the box
  
ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_small_centroid, mapping = aes (color = TDV))
  
  # Create a density that is the *exact* opposite of the natural density created by distances
  # We are basically forcing a near-perfect density homogeneity
  mu <-  mean(transition.distance.vector) / transition.distance.vector
  
  # Transition matrix with distances factored with mu 
  transition.distance.mu <- matrix(0, 400, 400) 
  
  for (i in 1: nrow(sf_grid_small)) { # For each row ("from")
    for (j in 1:nrow(sf_grid_small)) { # For each column ("to")
      transition.distance.mu[i,j] <- transition.distance[i,j] * mu[j] # Multiply distance with mu
    }
  }
  
  transition.distance.mu[1:20, 1:20]
  
  colSums(transition.distance.mu) # Now you can see that each cell has equal probability of being "moved to".
  
  # Create a dataframe for sampling probability // Northing will be sampling probability
  sampling.prob <- data.frame (ID.small = c(1:400), 
                               sampling.prob = rep (seq(from = 1, to = 10, by = 0.4736842), each = 20)) # Ridiculous by = n, but like that I have an exact 10-fold increase between first row and last row
  
  ###---Breeding---###
  
  ###---Year 0 breeding---###
  
  mothers <- which(init.pop.sf$sex=='F' & init.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
  
  YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
  
  # Loop through all the available mothers
  for(j in 1:length(mothers)){
    
    num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    
    # Record the mother
    mother <- as_Spatial(init.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
    
    # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
    fathers <- which(init.pop.sf$sex=='M' & init.pop.sf$age>=repro.age & init.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
    
    if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
    
    # Need to make a combination of number of available fathers and number of mates this year for this female
    ifelse(length(fathers) == 1,  # If only one father available...
           {num.mates.x <- 1; fathers.ID <- init.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
           ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
    
    num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
    num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
    
    if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
    
    # If there are two fathers, need to distribute the cubs among them
    ifelse (num.offspring %% 2 == 0, # If even number of cubs....
            num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
            num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
    
    for (h in 1:num.mates.x){
      
      ifelse (num.mates.x == 1, 
              num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
              ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                     num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                     num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
      
      if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
      
      indv.name <- NA # Create a place holder for the random name given to each offspring
      age <- 0 # Assign a 0 age to each offspring born this year
      birth.year <- 0 # Note that these pups were born in the 0 year of the simulation
      sex <- NA # Create a place holder for the sexes
      long <- NA
      lat <- NA
      dispersal <- NA
      
      
      inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                   "mother" = mother[[1]], # record the mother's name
                                   "father" = fathers.ID[[1]][[h]], # record the father's name,
                                   sex, long, lat, dispersal) 
      
      inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
      inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
      inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
      inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
      inter.df$dispersal <- 0 # Newborns have not dispersed yet
      
      baby.names <- vector() # Create a blank vector for random baby names
      
      for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
        name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
        baby.names <- c(baby.names,name1)
      } # End loop over name
      
      inter.df$indv.name <- baby.names # add the baby names to the data frame
      
      YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
      
    } # End loop over mates
    
  } # End loop over mothers
  
  YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
  
  # Assign each YOY to the cell they are into
  YOY.df.sf <- YOY.df.sf %>% 
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  
  # Rename columns for convenience
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'  
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
  
  
  loopy.pop.sf <- NULL
  loopy.pop.sf <- rbind(init.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
  
  # At end of year 0
  Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  
  # MORTALITY (two steps)
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
  loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                  ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                         ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
  
  # Density-dependent yoy survival
  next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
  diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
  ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
  yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
  loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
  
  # Lethal sampling of individuals. 
  
  # Size of this year's adult & juvenile population (before survival draw).
  number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
  
  # Add the northing (sampling probability) to population dataframe
  loopy.pop.sf$sampling.prob <- sampling.prob[as.numeric(loopy.pop.sf$ID.small),"sampling.prob"]
  
  sample.pop <- loopy.pop.sf %>%
    subset (survival == "M") %>% # Only harvest individuals
    subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
    slice_sample(n = round((harvest.size*number.adults)/100), weight_by = sampling.prob, replace = FALSE) # Harvest n% of this year's total adult population
  
  # Remove northing (I could also just leave the northing to the entire dataset every time..)
  loopy.pop.sf <- subset (loopy.pop.sf, select = - sampling.prob)
  
  # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
  loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                   "H", # Then write him as harvested
                                   loopy.pop.sf$survival) # Otherwise keep S or M as it was before
  
  # Dispersal and grid assignment 
  # Since all my initial pop has dispersed, and none of the YOY for year 0 disperse, there is no dispersal event for the first year
  
  # At end of year 0
  print(paste("year 0 ", "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
  
  # Note: not all males or females are fathers or mothers
  # Mothers and fathers may have contributed to reproduction this year and die afterwards (since survival is estimated at a later stage)
  
  # Visualize the YOY
  # Don't worry if it looks like there are more males than females. The sibling dots overlap, and males are "on top" of females, but sex ratio is okay.
  # ggplot() +
  #   geom_sf(data = sf_box, fill = 'white', lwd = 0.05) +
  #   geom_sf(data = sf_grid_small, fill = 'transparent', lwd = 0.3) +
  #   geom_sf(data = sample.pop, mapping = aes (color = sex), size = 0.7) #+
  #coord_sf(datum = st_crs(32616)) # To change axis to UTM (meters)
  
  ###---Loop for all other breeding years---###
  
  pop.size <- data.frame()
  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year
  sample.list <- list() # Make list to store dataframe of sampled individual for each year
  
  for(v in 1:(num.years)){ #loop through all of the years in the simulation
    
    v.pop.sf <- loopy.pop.sf[loopy.pop.sf$survival =="S",] # Bring in the data from the previous iteration, but only include those that survive
    v.pop.sf <- subset (v.pop.sf, select = - survival) # Remove the survival column
    v.pop.sf$age <- v.pop.sf$age+1 # Increase each individuals age by one for the new year
    
    mothers <- which(v.pop.sf$sex=='F' & v.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
    
    YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
    
    # Loop through all the available mothers
    for(j in 1:length(mothers)){
      
      num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
      inter.df <- data.frame() # Create a data frame for the offspring born from each mother
      
      # Record the mother
      mother <- as_Spatial(v.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
      
      # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
      fathers <- which(v.pop.sf$sex=='M' & v.pop.sf$age>=repro.age & v.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
      
      if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
      
      # Need to make a combination of number of available fathers and number of mates this year for this female
      ifelse(length(fathers) == 1,  # If only one father available...
             {num.mates.x <- 1; fathers.ID <- v.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
             ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
      
      num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
      num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
      
      if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
      
      # If there are two fathers, need to distribute the cubs among them
      ifelse (num.offspring %% 2 == 0, # If even number of cubs....
              num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
              num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
      
      for (h in 1:num.mates.x){
        
        ifelse (num.mates.x == 1, 
                num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
                ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                       num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                       num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
        
        if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
        
        indv.name <- NA # Create a place holder for the random name given to each offspring
        age <- 0 # Assign a 0 age to each offspring born this year
        birth.year <- v # Note that these pups were born in the 0 year of the simulation
        sex <- NA # Create a place holder for the sexes
        long <- NA
        lat <- NA
        dispersal <- NA
        
        
        inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                     "mother" = mother[[1]], # record the mother's name
                                     "father" = fathers.ID[[1]][[h]], # record the father's name,
                                     sex, long, lat, dispersal) 
        
        inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
        inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
        inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
        inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
        inter.df$dispersal <- 0 # Newborns have not dispersed yet
        
        baby.names <- vector() # Create a blank vector for random baby names
        
        for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
          name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
          baby.names <- c(baby.names,name1)
        } # End loop over name
        
        inter.df$indv.name <- baby.names # add the baby names to the data frame
        
        YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
        
      } # End loop over mates
      
    } # End loop over mothers
    
    YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
    
    # Assign each YOY to the cell they are into
    YOY.df.sf <- YOY.df.sf %>% 
      st_join(sf_grid_small, join = st_intersects)  %>% # Spatial join of each point to the grid ID. Adds a new column.
      st_join(sf_grid_large, join = st_intersects)
    
    
    # Rename columns for convenience
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
    
    loopy.pop.sf <- NULL
    loopy.pop.sf <- rbind(v.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
    
    # At end of year 0
    Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    
    # DISPERSAL
    
    # Let's fix dispersal at 1yo. 
    dispersers.pop <- data.frame() # Create a place holder for this year's dispersers
    
    dispersers.pop <- loopy.pop.sf %>%
      subset (dispersal == 0) %>% # Only one dispersal events, so remove those that have already dispersed (only used if eg only 50% of 1yo dispersed, and the rest at 2yo)
      subset (age == 1) %>% # All 1 year olds disperse
      st_drop_geometry() %>% # Remove geometry for now
      subset (select = - ID.large) # Remove large cell ID (it might change after dispersal)
    
    for (d in 1:nrow(dispersers.pop)) {
      
      dest.vec <- transition.distance.mu[as.numeric(dispersers.pop[d, "ID.small"]),] # Get the cell # of the original location ("from)
      
      dispersers.pop[d, "ID.small"] <- sample(1:nrow(sf_grid_small), size = 1, prob = dest.vec) # Select final destination ("to") based on distance x mu (probability of same point = 0)
      
      # Retrieve geometry from small grid ID (based on centroids), transformed into long and lat columns
      dispersers.pop[d, "long"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "X"]
      dispersers.pop[d, "lat"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "Y"]
    }
    
    dispersers.pop.sf <- st_as_sf(dispersers.pop, coords = c("long", "lat"), crs = 4326) %>%  # Create the spatial feature for dispersers
      st_join(sf_grid_large, join = st_intersects) # Spatial join of each point to the large grid ID. Adds a new column.
    
    # Rename column for convenience
    names(dispersers.pop.sf)[names(dispersers.pop.sf) == 'ID'] <- 'ID.large'
    
    # Assign dispersal as being done
    dispersers.pop.sf$dispersal <- 1
    
    # DIAGNOSTIC: let's connect the "start" and "end" of dispersal. Gives an idea of how far individuals move around.
    # 
    # dispersers.pop.origin <- loopy.pop.sf %>% # This one is just for diagnostic purposes
    #   subset (dispersal == 0) %>%
    #   subset (age == 1)
    # 
    # start.long <- st_coordinates(dispersers.pop.origin)[,1]
    # start.lat <- st_coordinates(dispersers.pop.origin)[,2]
    # 
    # end.long <- st_coordinates(dispersers.pop.sf)[,1]
    # end.lat <- st_coordinates(dispersers.pop.sf)[,2]
    # 
    # # plot(c(start.long, end.long), c(start.lat, end.lat), col = rep(c("red", "blue"), each = length(start.long)),
    # #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # # arrows(start.long, start.lat, end.long, end.lat)
    # #
    # # # If it's too messy, you can subset the dataframes
    # 
    # start.long.sub <- start.long[1:(length(start.long)/5)]
    # start.lat.sub <- start.lat[1:(length(start.lat)/5)]
    # 
    # end.long.sub <- end.long[1:(length(end.long)/5)]
    # end.lat.sub <- end.lat[1:(length(end.lat)/5)]
    # 
    # plot(c(start.long.sub, end.long.sub), c(start.lat.sub, end.lat.sub), col = rep(c("red", "blue"), each = length(start.long.sub)),
    #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # arrows(start.long.sub, start.lat.sub, end.long.sub, end.lat.sub)
    # 
    # Remove individuals that are in dispersers.pop from the loopy.pop, and add them again with the new dataset
    loopy.pop.sf <- loopy.pop.sf[!loopy.pop.sf$indv.name %in% dispersers.pop.sf$indv.name,] # Remove
    loopy.pop.sf <- rbind.data.frame(loopy.pop.sf, dispersers.pop.sf) # Add
    
    # MORTALITY (two steps)
    
    #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
    loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                    ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                           ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
    
    # Density-dependent yoy survival
    next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
    diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
    ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
    yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
    loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
    
    # Lethal sampling of individuals. 
    
    # Size of this year's adult & juvenile population (before survival draw).
    number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
    
    # Add the northing (sampling probability) to population dataframe
    loopy.pop.sf$sampling.prob <- sampling.prob[as.numeric(loopy.pop.sf$ID.small),"sampling.prob"]
    
    sample.pop <- loopy.pop.sf %>%
      subset (survival == "M") %>% # Only harvest individuals
      subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
      slice_sample(n = round((harvest.size*number.adults)/100), weight_by = sampling.prob, replace = FALSE) # Harvest n% of this year's total adult population
    
    # Remove northing
    loopy.pop.sf <- subset (loopy.pop.sf, select = - sampling.prob)
    
    # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
    loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                     "H", # Then write him as harvested
                                     loopy.pop.sf$survival) # Otherwise keep S or M as it was before
    
    
    # Return samples to a data frame, with geometry transformed into long and lat columns
    sample.pop <- sample.pop %>%
      mutate(long = st_coordinates(.)[,"X"], #Re-extract long and lat coordinates in separate columns
             lat = st_coordinates(.)[,"Y"]) %>%
      st_set_geometry(NULL) # Remove geometry
    
    sample.list[[v]] <- sample.pop
    loopy.list[[v]] <- loopy.pop.sf # Save the current year's population data as a list element, where the index corresponds to the year
    
    # Print the simulation year and the population size in the R console so they can be observed
    print(paste("year", v, "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df.sf), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
    
    ### Count survivors for each year v
    pop.size.vec <- cbind.data.frame(year=v,
                                     population_size_preharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]), # Minus young of the year
                                     population_size_postharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1 & loopy.pop.sf$survival=="S",]), # Minus young of the year
                                     Nm_preharvest = Nm, # Includes males that may have reproduced this year but were subsequently killed (harvest or natural)
                                     Nf_preharvest = Nf) # Includes females that may have reproduced this year but were subsequently killed (harvest or natural)
    
    pop.size <- rbind(pop.size, pop.size.vec)
    
  } # end loop over sim years
  
  # Label the list elements with the year
  names(loopy.list) <- paste0("year.end.pop.", seq(1:num.years))
  
  return(invisible(list(loopy.list, sample.list, pop.size)))
  
}


HPP.restricted <- function(){
  
  # This is an attempt to "force" population to be homogeneous (i.e., force some dispersal on the outside of the box)
  
  ###---Set up your initial population---###
  
  init.pop <- data.frame() # Create a blank data frame which will become the initial population
  
  age <- NULL # Create an empty vector for age
  for(y in 1:length(Nages)) age <- c(age, rep(y, Nages[y])) # Elongate age into population size and fill each row according to ratios in Nages
  
  # Loop that creates the below data for each individual in a population the size of "init.pop.size"
  for(i in 1:init.pop.size){
    
    indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
    mother <- "xxxxx" # The individuals in the initial population do not have known mothers
    father <- "xxxxx" # The individuals in the initial population do not have known fathers
    birth.year <- -1 # this is a place holder for individuals born within the simulation
    sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the proportions set in the parameter section
    dispersal <- 1 # Let's assume they all dispersed, no matter the age (only 1 dispersal event in life)
    init.vec <- cbind.data.frame(indv.name, birth.year, age[i], mother, father, sex, dispersal) # Create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
    init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
  }
  
  names(init.pop)[3] <- "age" # Change the name to remove the [i] from the loop
  
  # Have a look at your population
  
  head(init.pop)
  summary(as.factor(init.pop$age)) # Should match Nages
  summary(as.factor(init.pop$sex)) # Should be 50/50
  
  # Create a box (100x100km)
  box <- matrix(c(0, 0, # Create the bounding points
                  0, 1,
                  1, 1,
                  1, 0,
                  0, 0),
                ncol=2, byrow=TRUE)
  
  sf_box <- st_polygon(list(box)) # Join the points to create a polygon
  sf_box <- st_sfc(sf_box, crs = 4326) # Transform polygon into a spatial feature
  
  # Then create two grids on box
  # One larger grid for reproduction
  sf_grid_large  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.1, 0.1)) %>% # Make it 10x10
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326) # Reverts the dataframe to a sf object (now including ID)
  
  # One smaller grid for dispersal and heterogenous point process (/density)
  sf_grid_small  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.05, 0.05)) %>% # Make it 50x50
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326)
  
  # Get centroid for smaller grid
  sf_grid_small_centroid <- st_centroid(sf_grid_small)
  
  # Create random points that will populate the small grid (centroid of each cell)
  points_init <- st_sample(sf_grid_small_centroid, size = nrow(init.pop), replace = TRUE)
  
  # Assign those points to the init.pop dataframe
  init.pop$long <- st_coordinates(points_init)[,1]
  init.pop$lat <- st_coordinates(points_init)[,2]
  
  # Assign each ind to the cell they are into
  init.pop.sf <- st_as_sf (init.pop, coords = c("long", "lat"), crs = 4326) %>% # Transform your dataframe as a sf
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  # Rename columns for convenience
  names(init.pop.sf)[names(init.pop.sf) == 'ID.x'] <- 'ID.small'
  names(init.pop.sf)[names(init.pop.sf) == 'ID.y'] <- 'ID.large'
  
  # Visualize what you just did
  
  sf_grid_small_centroid.coord <- as.data.frame(sf::st_coordinates(sf_grid_small_centroid))
  sf_grid_small_centroid.coord$ID <- sf_grid_small_centroid$ID
  
  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_large, lwd = 0.3) +
    geom_sf(data = sf_grid_small, lwd = 0.3) +
    #geom_sf(data = init.pop.sf, mapping = aes (color = sex), size = 0.7) +
    # geom_sf(data = sf_grid_small_centroid, size = 0.7) +
    geom_text(data = sf_grid_small_centroid.coord, aes(X, Y, label = ID), colour = "black")
  
  # Create matrix and fill with distances
  distance.matrix <- st_distance(x = sf_grid_small_centroid, y = sf_grid_small_centroid) # units = m
  distance.matrix[1:10, 1:10] # ~5.5km between adjacent grids
  
  # Transition matrix based on distances only (follows dgamma with shape = 2 and rate = 1.5e-4)
  # Note: might need to do one matrix for males and one matrix for females if dispersal distances are different
  transition.distance = matrix(dgamma(as.numeric(distance.matrix), 2, 1.5e-4), 400, 400)
  transition.distance[1:10, 1:10]
  
  # Divide by a normalizing constant so that sum(transition.distance) == 1
  transition.distance <- transition.distance / sum(transition.distance)
  
  # Calculate SUM probability of ind moving to each cell based on distance only
  transition.distance.vector <- colSums(transition.distance) 
  sf_grid_small_centroid$TDV <- transition.distance.vector # Using the plot below, you can see that there will naturally be more points in the center of the box
  
  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_small_centroid, mapping = aes (color = TDV))
  
  # Create a density that is the *exact* opposite of the natural density created by distances
  # We are basically forcing a near-perfect density homogeneity
  mu <-  mean(transition.distance.vector) / transition.distance.vector
  
  # Transition matrix with distances factored with mu 
  transition.distance.mu <- matrix(0, 400, 400) 
  
  for (i in 1: nrow(sf_grid_small)) { # For each row ("from")
    for (j in 1:nrow(sf_grid_small)) { # For each column ("to")
      transition.distance.mu[i,j] <- transition.distance[i,j] * mu[j] # Multiply distance with mu
    }
  }
  
  transition.distance.mu[1:20, 1:20]
  
  colSums(transition.distance.mu) # Now you can see that each cell has equal probability of being "moved to".
  
  # Create a dataframe for sampling probability // Northing will be sampling probability
  sampling.prob <- data.frame (ID.small = c(1:400), 
                               sampling.prob = c(rep(0, times = 200),
                                                 rep(seq(from = 1, to = 10, by = 1), each = 20)))
  ###---Breeding---###
  
  ###---Year 0 breeding---###
  
  mothers <- which(init.pop.sf$sex=='F' & init.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
  
  YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
  
  # Loop through all the available mothers
  for(j in 1:length(mothers)){
    
    num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    
    # Record the mother
    mother <- as_Spatial(init.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
    
    # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
    fathers <- which(init.pop.sf$sex=='M' & init.pop.sf$age>=repro.age & init.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
    
    if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
    
    # Need to make a combination of number of available fathers and number of mates this year for this female
    ifelse(length(fathers) == 1,  # If only one father available...
           {num.mates.x <- 1; fathers.ID <- init.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
           ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
    
    num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
    num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
    
    if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
    
    # If there are two fathers, need to distribute the cubs among them
    ifelse (num.offspring %% 2 == 0, # If even number of cubs....
            num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
            num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
    
    for (h in 1:num.mates.x){
      
      ifelse (num.mates.x == 1, 
              num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
              ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                     num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                     num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
      
      if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
      
      indv.name <- NA # Create a place holder for the random name given to each offspring
      age <- 0 # Assign a 0 age to each offspring born this year
      birth.year <- 0 # Note that these pups were born in the 0 year of the simulation
      sex <- NA # Create a place holder for the sexes
      long <- NA
      lat <- NA
      dispersal <- NA
      
      
      inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                   "mother" = mother[[1]], # record the mother's name
                                   "father" = fathers.ID[[1]][[h]], # record the father's name,
                                   sex, long, lat, dispersal) 
      
      inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
      inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
      inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
      inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
      inter.df$dispersal <- 0 # Newborns have not dispersed yet
      
      baby.names <- vector() # Create a blank vector for random baby names
      
      for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
        name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
        baby.names <- c(baby.names,name1)
      } # End loop over name
      
      inter.df$indv.name <- baby.names # add the baby names to the data frame
      
      YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
      
    } # End loop over mates
    
  } # End loop over mothers
  
  YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
  
  # Assign each YOY to the cell they are into
  YOY.df.sf <- YOY.df.sf %>% 
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  
  # Rename columns for convenience
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'  
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
  
  
  loopy.pop.sf <- NULL
  loopy.pop.sf <- rbind(init.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
  
  # At end of year 0
  Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  
  # MORTALITY (two steps)
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
  loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                  ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                         ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
  
  # Density-dependent yoy survival
  next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
  diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
  ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
  yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
  loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
  
  # Lethal sampling of individuals. 
  
  # Size of this year's adult & juvenile population (before survival draw).
  number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
  
  # Add the northing (sampling probability) to population dataframe
  loopy.pop.sf$sampling.prob <- sampling.prob[as.numeric(loopy.pop.sf$ID.small),"sampling.prob"]
  
  sample.pop <- loopy.pop.sf %>%
    subset (survival == "M") %>% # Only harvest individuals
    subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
    slice_sample(n = round((harvest.size*number.adults)/100), weight_by = sampling.prob, replace = FALSE) # Harvest n% of this year's total adult population
  
  # Remove northing (I could also just leave the northing to the entire dataset every time..)
  loopy.pop.sf <- subset (loopy.pop.sf, select = - sampling.prob)
  
  # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
  loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                   "H", # Then write him as harvested
                                   loopy.pop.sf$survival) # Otherwise keep S or M as it was before
  
  # Dispersal and grid assignment 
  # Since all my initial pop has dispersed, and none of the YOY for year 0 disperse, there is no dispersal event for the first year
  
  # At end of year 0
  print(paste("year 0 ", "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
  
  # Note: not all males or females are fathers or mothers
  # Mothers and fathers may have contributed to reproduction this year and die afterwards (since survival is estimated at a later stage)
  
  # Visualize the YOY
  # Don't worry if it looks like there are more males than females. The sibling dots overlap, and males are "on top" of females, but sex ratio is okay.
  # ggplot() +
  #   geom_sf(data = sf_box, fill = 'white', lwd = 0.05) +
  #   geom_sf(data = sf_grid_small, fill = 'transparent', lwd = 0.3) +
  #   geom_sf(data = sample.pop, mapping = aes (color = sex), size = 0.7) #+
  #coord_sf(datum = st_crs(32616)) # To change axis to UTM (meters)
  
  ###---Loop for all other breeding years---###
  
  pop.size <- data.frame()
  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year
  sample.list <- list() # Make list to store dataframe of sampled individual for each year
  
  for(v in 1:(num.years)){ #loop through all of the years in the simulation
    
    v.pop.sf <- loopy.pop.sf[loopy.pop.sf$survival =="S",] # Bring in the data from the previous iteration, but only include those that survive
    v.pop.sf <- subset (v.pop.sf, select = - survival) # Remove the survival column
    v.pop.sf$age <- v.pop.sf$age+1 # Increase each individuals age by one for the new year
    
    mothers <- which(v.pop.sf$sex=='F' & v.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
    
    YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
    
    # Loop through all the available mothers
    for(j in 1:length(mothers)){
      
      num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
      inter.df <- data.frame() # Create a data frame for the offspring born from each mother
      
      # Record the mother
      mother <- as_Spatial(v.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
      
      # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
      fathers <- which(v.pop.sf$sex=='M' & v.pop.sf$age>=repro.age & v.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
      
      if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
      
      # Need to make a combination of number of available fathers and number of mates this year for this female
      ifelse(length(fathers) == 1,  # If only one father available...
             {num.mates.x <- 1; fathers.ID <- v.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
             ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
      
      num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
      num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
      
      if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
      
      # If there are two fathers, need to distribute the cubs among them
      ifelse (num.offspring %% 2 == 0, # If even number of cubs....
              num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
              num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
      
      for (h in 1:num.mates.x){
        
        ifelse (num.mates.x == 1, 
                num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
                ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                       num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                       num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
        
        if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
        
        indv.name <- NA # Create a place holder for the random name given to each offspring
        age <- 0 # Assign a 0 age to each offspring born this year
        birth.year <- v # Note that these pups were born in the 0 year of the simulation
        sex <- NA # Create a place holder for the sexes
        long <- NA
        lat <- NA
        dispersal <- NA
        
        
        inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                     "mother" = mother[[1]], # record the mother's name
                                     "father" = fathers.ID[[1]][[h]], # record the father's name,
                                     sex, long, lat, dispersal) 
        
        inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
        inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
        inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
        inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
        inter.df$dispersal <- 0 # Newborns have not dispersed yet
        
        baby.names <- vector() # Create a blank vector for random baby names
        
        for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
          name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
          baby.names <- c(baby.names,name1)
        } # End loop over name
        
        inter.df$indv.name <- baby.names # add the baby names to the data frame
        
        YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
        
      } # End loop over mates
      
    } # End loop over mothers
    
    YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
    
    # Assign each YOY to the cell they are into
    YOY.df.sf <- YOY.df.sf %>% 
      st_join(sf_grid_small, join = st_intersects)  %>% # Spatial join of each point to the grid ID. Adds a new column.
      st_join(sf_grid_large, join = st_intersects)
    
    
    # Rename columns for convenience
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
    
    loopy.pop.sf <- NULL
    loopy.pop.sf <- rbind(v.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
    
    # At end of year 0
    Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    
    # DISPERSAL
    
    # Let's fix dispersal at 1yo. 
    dispersers.pop <- data.frame() # Create a place holder for this year's dispersers
    
    dispersers.pop <- loopy.pop.sf %>%
      subset (dispersal == 0) %>% # Only one dispersal events, so remove those that have already dispersed (only used if eg only 50% of 1yo dispersed, and the rest at 2yo)
      subset (age == 1) %>% # All 1 year olds disperse
      st_drop_geometry() %>% # Remove geometry for now
      subset (select = - ID.large) # Remove large cell ID (it might change after dispersal)
    
    for (d in 1:nrow(dispersers.pop)) {
      
      dest.vec <- transition.distance.mu[as.numeric(dispersers.pop[d, "ID.small"]),] # Get the cell # of the original location ("from)
      
      dispersers.pop[d, "ID.small"] <- sample(1:nrow(sf_grid_small), size = 1, prob = dest.vec) # Select final destination ("to") based on distance x mu (probability of same point = 0)
      
      # Retrieve geometry from small grid ID (based on centroids), transformed into long and lat columns
      dispersers.pop[d, "long"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "X"]
      dispersers.pop[d, "lat"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "Y"]
    }
    
    dispersers.pop.sf <- st_as_sf(dispersers.pop, coords = c("long", "lat"), crs = 4326) %>%  # Create the spatial feature for dispersers
      st_join(sf_grid_large, join = st_intersects) # Spatial join of each point to the large grid ID. Adds a new column.
    
    # Rename column for convenience
    names(dispersers.pop.sf)[names(dispersers.pop.sf) == 'ID'] <- 'ID.large'
    
    # Assign dispersal as being done
    dispersers.pop.sf$dispersal <- 1
    
    # DIAGNOSTIC: let's connect the "start" and "end" of dispersal. Gives an idea of how far individuals move around.
    #
    # dispersers.pop.origin <- loopy.pop.sf %>% # This one is just for diagnostic purposes
    #   subset (dispersal == 0) %>%
    #   subset (age == 1)
    # 
    # start.long <- st_coordinates(dispersers.pop.origin)[,1]
    # start.lat <- st_coordinates(dispersers.pop.origin)[,2]
    # 
    # end.long <- st_coordinates(dispersers.pop.sf)[,1]
    # end.lat <- st_coordinates(dispersers.pop.sf)[,2]
    # # 
    # # # plot(c(start.long, end.long), c(start.lat, end.lat), col = rep(c("red", "blue"), each = length(start.long)),
    # # #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # # # arrows(start.long, start.lat, end.long, end.lat)
    # # 
    # # # If it's too messy, you can subset the dataframes
    # # 
    # start.long.sub <- start.long[1:(length(start.long)/5)]
    # start.lat.sub <- start.lat[1:(length(start.lat)/5)]
    # 
    # end.long.sub <- end.long[1:(length(end.long)/5)]
    # end.lat.sub <- end.lat[1:(length(end.lat)/5)]
    # 
    # plot(c(start.long.sub, end.long.sub), c(start.lat.sub, end.lat.sub), col = rep(c("red", "blue"), each = length(start.long.sub)),
    #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # arrows(start.long.sub, start.lat.sub, end.long.sub, end.lat.sub)

    # Remove individuals that are in dispersers.pop from the loopy.pop, and add them again with the new dataset
    loopy.pop.sf <- loopy.pop.sf[!loopy.pop.sf$indv.name %in% dispersers.pop.sf$indv.name,] # Remove
    loopy.pop.sf <- rbind.data.frame(loopy.pop.sf, dispersers.pop.sf) # Add
    
    # MORTALITY (two steps)
    
    #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
    loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                    ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                           ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
    
    # Density-dependent yoy survival
    next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
    diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
    ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
    yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
    loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
    
    # Lethal sampling of individuals. 
    
    # Size of this year's adult & juvenile population (before survival draw).
    number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
    
    # Add the northing (sampling probability) to population dataframe
    loopy.pop.sf$sampling.prob <- sampling.prob[as.numeric(loopy.pop.sf$ID.small),"sampling.prob"]
    
    sample.pop <- loopy.pop.sf %>%
      subset (survival == "M") %>% # Only harvest individuals
      subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
      slice_sample(n = round((harvest.size*number.adults)/100), weight_by = sampling.prob, replace = FALSE) # Harvest n% of this year's total adult population
    
    # Remove northing
    loopy.pop.sf <- subset (loopy.pop.sf, select = - sampling.prob)
    
    # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
    loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                     "H", # Then write him as harvested
                                     loopy.pop.sf$survival) # Otherwise keep S or M as it was before
    
    
    # Return samples to a data frame, with geometry transformed into long and lat columns
    sample.pop <- sample.pop %>%
      mutate(long = st_coordinates(.)[,"X"], #Re-extract long and lat coordinates in separate columns
             lat = st_coordinates(.)[,"Y"]) %>%
      st_set_geometry(NULL) # Remove geometry
    
    sample.list[[v]] <- sample.pop
    loopy.list[[v]] <- loopy.pop.sf # Save the current year's population data as a list element, where the index corresponds to the year
    
    # Print the simulation year and the population size in the R console so they can be observed
    print(paste("year", v, "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df.sf), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
    
    ### Count survivors for each year v
    pop.size.vec <- cbind.data.frame(year=v,
                                     population_size_preharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]), # Minus young of the year
                                     population_size_postharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1 & loopy.pop.sf$survival=="S",]), # Minus young of the year
                                     Nm_preharvest = Nm, # Includes males that may have reproduced this year but were subsequently killed (harvest or natural)
                                     Nf_preharvest = Nf) # Includes females that may have reproduced this year but were subsequently killed (harvest or natural)
    
    pop.size <- rbind(pop.size, pop.size.vec)
    
  } # end loop over sim years
  
  # Label the list elements with the year
  names(loopy.list) <- paste0("year.end.pop.", seq(1:num.years))
  
  return(invisible(list(loopy.list, sample.list, pop.size)))
  
}

IPP.uniform <- function(){
  
  # Most basic population every (aka "baseline scenario", with equal and constant reproductive output)
  # Simple age at maturity, maximum age etc...
  
  ###---Set up your initial population---###
  
  init.pop <- data.frame() # Create a blank data frame which will become the initial population
  
  age <- NULL # Create an empty vector for age
  for(y in 1:length(Nages)) age <- c(age, rep(y, Nages[y])) # Elongate age into population size and fill each row according to ratios in Nages
  
  # Loop that creates the below data for each individual in a population the size of "init.pop.size"
  for(i in 1:init.pop.size){
    
    indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
    mother <- "xxxxx" # The individuals in the initial population do not have known mothers
    father <- "xxxxx" # The individuals in the initial population do not have known fathers
    birth.year <- -1 # this is a place holder for individuals born within the simulation
    sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the proportions set in the parameter section
    dispersal <- 1 # Let's assume they all dispersed, no matter the age (only 1 dispersal event in life)
    init.vec <- cbind.data.frame(indv.name, birth.year, age[i], mother, father, sex, dispersal) # Create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
    init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
  }
  
  names(init.pop)[3] <- "age" # Change the name to remove the [i] from the loop
  
  # Have a look at your population
  
  head(init.pop)
  summary(as.factor(init.pop$age)) # Should match Nages
  summary(as.factor(init.pop$sex)) # Should be 50/50
  
  # Create a box (100x100km)
  box <- matrix(c(0, 0, # Create the bounding points
                  0, 1,
                  1, 1,
                  1, 0,
                  0, 0),
                ncol=2, byrow=TRUE)
  
  sf_box <- st_polygon(list(box)) # Join the points to create a polygon
  sf_box <- st_sfc(sf_box, crs = 4326) # Transform polygon into a spatial feature
  
  # Then create two grids on box
  # One larger grid for reproduction
  sf_grid_large  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.1, 0.1)) %>% # Make it 10x10
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326) # Reverts the dataframe to a sf object (now including ID)
  
  # One smaller grid for dispersal and heterogenous point process (/density)
  sf_grid_small  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.05, 0.05)) %>% # Make it 50x50
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326)
  
  # Get centroid for smaller grid
  sf_grid_small_centroid <- st_centroid(sf_grid_small)
  
  # Create random points that will populate the small grid (centroid of each cell)
  points_init <- st_sample(sf_grid_small_centroid, size = nrow(init.pop), replace = TRUE)
  
  # Assign those points to the init.pop dataframe
  init.pop$long <- st_coordinates(points_init)[,1]
  init.pop$lat <- st_coordinates(points_init)[,2]
  
  # Assign each ind to the cell they are into
  init.pop.sf <- st_as_sf (init.pop, coords = c("long", "lat"), crs = 4326) %>% # Transform your dataframe as a sf
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  # Rename columns for convenience
  names(init.pop.sf)[names(init.pop.sf) == 'ID.x'] <- 'ID.small'
  names(init.pop.sf)[names(init.pop.sf) == 'ID.y'] <- 'ID.large'
  
  # Visualize what you just did
  
  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_large, lwd = 0.3) +
    geom_sf(data = sf_grid_small, lwd = 0.3) +
    geom_sf(data = init.pop.sf, mapping = aes (color = sex), size = 0.7)
  
  # Create matrix and fill with distances
  distance.matrix <- st_distance(x = sf_grid_small_centroid, y = sf_grid_small_centroid) # units = m
  distance.matrix[1:10, 1:10] # ~5.5km between adjacent grids
  
  # Transition matrix based on distances only (follows dgamma)
  # Note: might need to do one matrix for males and one matrix for females if dispersal distances are different
  transition.distance = matrix(dgamma(as.numeric(distance.matrix), 2, 1.5e-4), 400, 400)
  
  # Divide by a normalizing constant so that sum(transition.distance) == 1
  transition.distance <- transition.distance / sum(transition.distance)
  
  # Calculate SUM probability of ind moving to each cell based on distance only
  transition.distance.vector <- colSums(transition.distance)
  sf_grid_small_centroid$TDV <- transition.distance.vector # Using the plot below, you can see that there will naturally be more points in the center of the box

  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_small_centroid, mapping = aes (color = TDV))

  # Create a density that is the *exact* opposite of the natural density created by distances
  mu_opposite <-  mean(transition.distance.vector) / transition.distance.vector

  # Create the intensity parameter mu
  mu <- c (rep(seq(from = 10, to = 1, by = -1), each = 20),
           rep(seq(from = 1, to = 10, by = 1), each = 20))

  # Empty matrix
  transition.distance.mu <- matrix(0, 400, 400)

  for (i in 1: nrow(sf_grid_small)) { # For each row ("from")
    for (j in 1:nrow(sf_grid_small)) { # For each column ("to")
      transition.distance.mu[i,j] <- transition.distance[i,j] * mu[j] * mu_opposite[j] # Multiply distance with mu and the mu_opposite
    }
  }

  mu <- colSums(transition.distance.mu) # Now you can see that cells in the same row have equal probability of being "moved to", whilst having the IPP

  # Diagnostic
  max(mu) / min(mu) # We are not looking at max = 10 and min = 1 anymore because of the factor by my_opposite, BUT the ratio between max and min should still be 10

  
  
  
  ###---Breeding---###
  
  ###---Year 0 breeding---###
  
  mothers <- which(init.pop.sf$sex=='F' & init.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
  
  YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
  
  # Loop through all the available mothers
  for(j in 1:length(mothers)){
    
    num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    
    # Record the mother
    mother <- as_Spatial(init.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
    
    # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
    fathers <- which(init.pop.sf$sex=='M' & init.pop.sf$age>=repro.age & init.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
    
    if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
    
    # Need to make a combination of number of available fathers and number of mates this year for this female
    ifelse(length(fathers) == 1,  # If only one father available...
           {num.mates.x <- 1; fathers.ID <- init.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
           ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
    
    num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
    num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
    
    if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
    
    # If there are two fathers, need to distribute the cubs among them
    ifelse (num.offspring %% 2 == 0, # If even number of cubs....
            num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
            num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
    
    for (h in 1:num.mates.x){
      
      ifelse (num.mates.x == 1, 
              num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
              ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                     num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                     num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
      
      if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
      
      indv.name <- NA # Create a place holder for the random name given to each offspring
      age <- 0 # Assign a 0 age to each offspring born this year
      birth.year <- 0 # Note that these pups were born in the 0 year of the simulation
      sex <- NA # Create a place holder for the sexes
      long <- NA
      lat <- NA
      dispersal <- NA
      
      
      inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                   "mother" = mother[[1]], # record the mother's name
                                   "father" = fathers.ID[[1]][[h]], # record the father's name,
                                   sex, long, lat, dispersal) 
      
      inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
      inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
      inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
      inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
      inter.df$dispersal <- 0 # Newborns have not dispersed yet
      
      baby.names <- vector() # Create a blank vector for random baby names
      
      for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
        name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
        baby.names <- c(baby.names,name1)
      } # End loop over name
      
      inter.df$indv.name <- baby.names # add the baby names to the data frame
      
      YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
      
    } # End loop over mates
    
  } # End loop over mothers
  
  YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
  
  # Assign each YOY to the cell they are into
  YOY.df.sf <- YOY.df.sf %>% 
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  
  # Rename columns for convenience
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'  
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
  
  
  loopy.pop.sf <- NULL
  loopy.pop.sf <- rbind(init.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
  
  # At end of year 0
  Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  
  # MORTALITY (two steps)
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
  loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                  ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                         ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
  
  # Density-dependent yoy survival
  next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
  diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
  ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
  yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
  loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
  
  # Lethal sampling of individuals. 
  
  # Size of this year's adult & juvenile population (before survival draw).
  number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
  
  sample.pop <- loopy.pop.sf %>%
    subset (survival == "M") %>% # Only harvest individuals
    subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
    slice_sample(n = round((harvest.size*number.adults)/100), replace = FALSE) # Harvest n% of this year's total adult population
  
  # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
  loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                   "H", # Then write him as harvested
                                   loopy.pop.sf$survival) # Otherwise keep S or M as it was before
  
  # Dispersal and grid assignment 
  # Since all my initial pop has dispersed, and none of the YOY for year 0 disperse, there is no dispersal event for the first year
  
  # At end of year 0
  print(paste("year 0 ", "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
  
  # Note: not all males or females are fathers or mothers
  # Mothers and fathers may have contributed to reproduction this year and die afterwards (since survival is estimated at a later stage)
  
  # Visualize the YOY
  # Don't worry if it looks like there are more males than females. The sibling dots overlap, and males are "on top" of females, but sex ratio is okay.
  # ggplot() +
  #   geom_sf(data = sf_box, fill = 'white', lwd = 0.05) +
  #   geom_sf(data = sf_grid, fill = 'transparent', lwd = 0.3) +
  #   geom_sf(data = YOY.df.sf, mapping = aes (color = sex), size = 0.7) #+
  #coord_sf(datum = st_crs(32616)) # To change axis to UTM (meters)
  
  ###---Loop for all other breeding years---###
  
  pop.size <- data.frame()
  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year
  sample.list <- list() # Make list to store dataframe of sampled individual for each year
  
  for(v in 1:(num.years)){ #loop through all of the years in the simulation
    
    v.pop.sf <- loopy.pop.sf[loopy.pop.sf$survival =="S",] # Bring in the data from the previous iteration, but only include those that survive
    v.pop.sf <- subset (v.pop.sf, select = - survival) # Remove the survival column
    v.pop.sf$age <- v.pop.sf$age+1 # Increase each individuals age by one for the new year
    
    mothers <- which(v.pop.sf$sex=='F' & v.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
    
    YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
    
    # Loop through all the available mothers
    for(j in 1:length(mothers)){
      
      num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
      inter.df <- data.frame() # Create a data frame for the offspring born from each mother
      
      # Record the mother
      mother <- as_Spatial(v.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
      
      # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
      fathers <- which(v.pop.sf$sex=='M' & v.pop.sf$age>=repro.age & v.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
      
      if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
      
      # Need to make a combination of number of available fathers and number of mates this year for this female
      ifelse(length(fathers) == 1,  # If only one father available...
             {num.mates.x <- 1; fathers.ID <- v.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
             ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
      
      num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
      num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
      
      if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
      
      # If there are two fathers, need to distribute the cubs among them
      ifelse (num.offspring %% 2 == 0, # If even number of cubs....
              num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
              num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
      
      for (h in 1:num.mates.x){
        
        ifelse (num.mates.x == 1, 
                num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
                ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                       num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                       num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
        
        if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
        
        indv.name <- NA # Create a place holder for the random name given to each offspring
        age <- 0 # Assign a 0 age to each offspring born this year
        birth.year <- v # Note that these pups were born in the 0 year of the simulation
        sex <- NA # Create a place holder for the sexes
        long <- NA
        lat <- NA
        dispersal <- NA
        
        
        inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                     "mother" = mother[[1]], # record the mother's name
                                     "father" = fathers.ID[[1]][[h]], # record the father's name,
                                     sex, long, lat, dispersal) 
        
        inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
        inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
        inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
        inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
        inter.df$dispersal <- 0 # Newborns have not dispersed yet
        
        baby.names <- vector() # Create a blank vector for random baby names
        
        for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
          name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
          baby.names <- c(baby.names,name1)
        } # End loop over name
        
        inter.df$indv.name <- baby.names # add the baby names to the data frame
        
        YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
        
      } # End loop over mates
      
    } # End loop over mothers
    
    YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
    
    # Assign each YOY to the cell they are into
    YOY.df.sf <- YOY.df.sf %>% 
      st_join(sf_grid_small, join = st_intersects)  %>% # Spatial join of each point to the grid ID. Adds a new column.
      st_join(sf_grid_large, join = st_intersects)
    
    
    # Rename columns for convenience
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
    
    loopy.pop.sf <- NULL
    loopy.pop.sf <- rbind(v.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
    
    # At end of year 0
    Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    
    # DISPERSAL
    
    # Let's fix dispersal at 1yo. 
    dispersers.pop <- data.frame() # Create a place holder for this year's dispersers
    
    dispersers.pop <- loopy.pop.sf %>%
      subset (dispersal == 0) %>% # Only one dispersal events, so remove those that have already dispersed (only used if eg only 50% of 1yo dispersed, and the rest at 2yo)
      subset (age == 1) %>% # All 1 year olds disperse
      st_drop_geometry() %>% # Remove geometry for now
      subset (select = - ID.large) # Remove large cell ID (it might change after dispersal)
    
    for (d in 1:nrow(dispersers.pop)) {
      
      dest.vec <- transition.distance.mu[as.numeric(dispersers.pop[d, "ID.small"]),] # Get the cell # of the original location ("from)
      
      dispersers.pop[d, "ID.small"] <- sample(1:nrow(sf_grid_small), size = 1, prob = dest.vec) # Select final destination ("to") based on distance x mu (probability of same point = 0)
      
      # Retrieve geometry from small grid ID (based on centroids), transformed into long and lat columns
      dispersers.pop[d, "long"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "X"]
      dispersers.pop[d, "lat"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "Y"]
    }
    
    dispersers.pop.sf <- st_as_sf(dispersers.pop, coords = c("long", "lat"), crs = 4326) %>%  # Create the spatial feature for dispersers
      st_join(sf_grid_large, join = st_intersects) # Spatial join of each point to the large grid ID. Adds a new column.
    
    # Rename column for convenience
    names(dispersers.pop.sf)[names(dispersers.pop.sf) == 'ID'] <- 'ID.large'
    
    # Assign dispersal as being done
    dispersers.pop.sf$dispersal <- 1
    
    # DIAGNOSTIC: let's connect the "start" and "end" of dispersal. Gives an idea of how far individuals move around.
    #
    # dispersers.pop.origin <- loopy.pop.sf %>% # This one is just for diagnostic purposes
    #   subset (dispersal == 0) %>%
    #   subset (age == 1)
    # 
    # start.long <- st_coordinates(dispersers.pop.origin)[,1]
    # start.lat <- st_coordinates(dispersers.pop.origin)[,2]
    # 
    # end.long <- st_coordinates(dispersers.pop.sf)[,1]
    # end.lat <- st_coordinates(dispersers.pop.sf)[,2]
    # 
    # # plot(c(start.long, end.long), c(start.lat, end.lat), col = rep(c("red", "blue"), each = length(start.long)),
    # #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # # arrows(start.long, start.lat, end.long, end.lat)
    # 
    # # If it's too messy, you can subset the dataframes
    # 
    # start.long.sub <- start.long[1:(length(start.long)/5)]
    # start.lat.sub <- start.lat[1:(length(start.lat)/5)]
    # 
    # end.long.sub <- end.long[1:(length(end.long)/5)]
    # end.lat.sub <- end.lat[1:(length(end.lat)/5)]
    # 
    # plot(c(start.long.sub, end.long.sub), c(start.lat.sub, end.lat.sub), col = rep(c("red", "blue"), each = length(start.long.sub)),
    #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # arrows(start.long.sub, start.lat.sub, end.long.sub, end.lat.sub)
    
    # Remove individuals that are in dispersers.pop from the loopy.pop, and add them again with the new dataset
    loopy.pop.sf <- loopy.pop.sf[!loopy.pop.sf$indv.name %in% dispersers.pop.sf$indv.name,] # Remove
    loopy.pop.sf <- rbind.data.frame(loopy.pop.sf, dispersers.pop.sf) # Add
    
    # MORTALITY (two steps)
    
    #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
    loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                    ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                           ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
    
    # Density-dependent yoy survival
    next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
    diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
    ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
    yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
    loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
    
    # Lethal sampling of individuals. 
    
    # Size of this year's adult & juvenile population (before survival draw).
    number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
    
    sample.pop <- loopy.pop.sf %>%
      subset (survival == "M") %>% # Only harvest individuals
      subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
      slice_sample(n = round((harvest.size*number.adults)/100), replace = FALSE) # Harvest n% of this year's total adult population
    
    # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
    loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                     "H", # Then write him as harvested
                                     loopy.pop.sf$survival) # Otherwise keep S or M as it was before
    
    
    # Return samples to a data frame, with geometry transformed into long and lat columns
    sample.pop <- sample.pop %>%
      mutate(long = st_coordinates(.)[,"X"], #Re-extract long and lat coordinates in separate columns
             lat = st_coordinates(.)[,"Y"]) %>%
      st_set_geometry(NULL) # Remove geometry
    
    sample.list[[v]] <- sample.pop
    loopy.list[[v]] <- loopy.pop.sf # Save the current year's population data as a list element, where the index corresponds to the year
    
    # Print the simulation year and the population size in the R console so they can be observed
    print(paste("year", v, "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df.sf), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
    
    ### Count survivors for each year v
    pop.size.vec <- cbind.data.frame(year=v,
                                     population_size_preharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]), # Minus young of the year
                                     population_size_postharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1 & loopy.pop.sf$survival=="S",]), # Minus young of the year
                                     Nm_preharvest = Nm, # Includes males that may have reproduced this year but were subsequently killed (harvest or natural)
                                     Nf_preharvest = Nf) # Includes females that may have reproduced this year but were subsequently killed (harvest or natural)
    
    pop.size <- rbind(pop.size, pop.size.vec)
    
  } # end loop over sim years
  
  # Label the list elements with the year
  names(loopy.list) <- paste0("year.end.pop.", seq(1:num.years))
  
  return(invisible(list(loopy.list, sample.list, pop.size)))
  
}

IPP.gradient <- function(){
  
  # Most basic population every (aka "baseline scenario", with equal and constant reproductive output)
  # Simple age at maturity, maximum age etc...
  
  ###---Set up your initial population---###
  
  init.pop <- data.frame() # Create a blank data frame which will become the initial population
  
  age <- NULL # Create an empty vector for age
  for(y in 1:length(Nages)) age <- c(age, rep(y, Nages[y])) # Elongate age into population size and fill each row according to ratios in Nages
  
  # Loop that creates the below data for each individual in a population the size of "init.pop.size"
  for(i in 1:init.pop.size){
    
    indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
    mother <- "xxxxx" # The individuals in the initial population do not have known mothers
    father <- "xxxxx" # The individuals in the initial population do not have known fathers
    birth.year <- -1 # this is a place holder for individuals born within the simulation
    sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the proportions set in the parameter section
    dispersal <- 1 # Let's assume they all dispersed, no matter the age (only 1 dispersal event in life)
    init.vec <- cbind.data.frame(indv.name, birth.year, age[i], mother, father, sex, dispersal) # Create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
    init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
  }
  
  names(init.pop)[3] <- "age" # Change the name to remove the [i] from the loop
  
  # Have a look at your population
  
  head(init.pop)
  summary(as.factor(init.pop$age)) # Should match Nages
  summary(as.factor(init.pop$sex)) # Should be 50/50
  
  # Create a box (100x100km)
  box <- matrix(c(0, 0, # Create the bounding points
                  0, 1,
                  1, 1,
                  1, 0,
                  0, 0),
                ncol=2, byrow=TRUE)
  
  sf_box <- st_polygon(list(box)) # Join the points to create a polygon
  sf_box <- st_sfc(sf_box, crs = 4326) # Transform polygon into a spatial feature
  
  # Then create two grids on box
  # One larger grid for reproduction
  sf_grid_large  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.1, 0.1)) %>% # Make it 10x10
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326) # Reverts the dataframe to a sf object (now including ID)
  
  # One smaller grid for dispersal and heterogenous point process (/density)
  sf_grid_small  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.05, 0.05)) %>% # Make it 50x50
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326)
  
  # Get centroid for smaller grid
  sf_grid_small_centroid <- st_centroid(sf_grid_small)
  
  # Create random points that will populate the small grid (centroid of each cell)
  points_init <- st_sample(sf_grid_small_centroid, size = nrow(init.pop), replace = TRUE)
  
  # Assign those points to the init.pop dataframe
  init.pop$long <- st_coordinates(points_init)[,1]
  init.pop$lat <- st_coordinates(points_init)[,2]
  
  # Assign each ind to the cell they are into
  init.pop.sf <- st_as_sf (init.pop, coords = c("long", "lat"), crs = 4326) %>% # Transform your dataframe as a sf
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  # Rename columns for convenience
  names(init.pop.sf)[names(init.pop.sf) == 'ID.x'] <- 'ID.small'
  names(init.pop.sf)[names(init.pop.sf) == 'ID.y'] <- 'ID.large'
  
  # Visualize what you just did
  
  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_large, lwd = 0.3) +
    geom_sf(data = sf_grid_small, lwd = 0.3) +
    geom_sf(data = init.pop.sf, mapping = aes (color = sex), size = 0.7)
  
  # Create matrix and fill with distances
  distance.matrix <- st_distance(x = sf_grid_small_centroid, y = sf_grid_small_centroid) # units = m
  distance.matrix[1:10, 1:10] # ~5.5km between adjacent grids
  
  # Transition matrix based on distances only (follows dgamma with shape = 4 and rate = 1.5e-4)
  # Note: might need to do one matrix for males and one matrix for females if dispersal distances are different
  transition.distance = matrix(dgamma(as.numeric(distance.matrix), 2, 1.5e-4), 400, 400)
  
  # Divide by a normalizing constant so that sum(transition.distance) == 1
  transition.distance <- transition.distance / sum(transition.distance)
  
  # Calculate SUM probability of ind moving to each cell based on distance only
  transition.distance.vector <- colSums(transition.distance) 
  sf_grid_small_centroid$TDV <- transition.distance.vector # Using the plot below, you can see that there will naturally be more points in the center of the box
  
  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_small_centroid, mapping = aes (color = TDV))
  
  # Create a density that is the *exact* opposite of the natural density created by distances
  mu_opposite <-  mean(transition.distance.vector) / transition.distance.vector
  
  # Create the intensity parameter mu
  mu <- c (rep(seq(from = 10, to = 1, by = -1), each = 20),
           rep(seq(from = 1, to = 10, by = 1), each = 20))
  
  # Empty matrix
  transition.distance.mu <- matrix(0, 400, 400) 
  
  for (i in 1: nrow(sf_grid_small)) { # For each row ("from")
    for (j in 1:nrow(sf_grid_small)) { # For each column ("to")
      transition.distance.mu[i,j] <- transition.distance[i,j] * mu[j] * mu_opposite[j] # Multiply distance with mu and the mu_opposite
    }
  }
  
  mu <- colSums(transition.distance.mu) # Now you can see that cells in the same row have equal probability of being "moved to", whilst having the IPP
  
  # Create a dataframe for sampling probability // Northing will be sampling probability
  sampling.prob <- data.frame (ID.small = c(1:400), 
                               sampling.prob = rep (seq(from = 1, to = 10, by = 0.4736842), each = 20)) # Ridiculous by = n, but like that I have an exact 10-fold increase between first row and last row
  
  ###---Breeding---###
  
  ###---Year 0 breeding---###
  
  mothers <- which(init.pop.sf$sex=='F' & init.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
  
  YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
  
  # Loop through all the available mothers
  for(j in 1:length(mothers)){
    
    num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    
    # Record the mother
    mother <- as_Spatial(init.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
    
    # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
    fathers <- which(init.pop.sf$sex=='M' & init.pop.sf$age>=repro.age & init.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
    
    if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
    
    # Need to make a combination of number of available fathers and number of mates this year for this female
    ifelse(length(fathers) == 1,  # If only one father available...
           {num.mates.x <- 1; fathers.ID <- init.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
           ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
    
    num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
    num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
    
    if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
    
    # If there are two fathers, need to distribute the cubs among them
    ifelse (num.offspring %% 2 == 0, # If even number of cubs....
            num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
            num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
    
    for (h in 1:num.mates.x){
      
      ifelse (num.mates.x == 1, 
              num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
              ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                     num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                     num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
      
      if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
      
      indv.name <- NA # Create a place holder for the random name given to each offspring
      age <- 0 # Assign a 0 age to each offspring born this year
      birth.year <- 0 # Note that these pups were born in the 0 year of the simulation
      sex <- NA # Create a place holder for the sexes
      long <- NA
      lat <- NA
      dispersal <- NA
      
      
      inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                   "mother" = mother[[1]], # record the mother's name
                                   "father" = fathers.ID[[1]][[h]], # record the father's name,
                                   sex, long, lat, dispersal) 
      
      inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
      inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
      inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
      inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
      inter.df$dispersal <- 0 # Newborns have not dispersed yet
      
      baby.names <- vector() # Create a blank vector for random baby names
      
      for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
        name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
        baby.names <- c(baby.names,name1)
      } # End loop over name
      
      inter.df$indv.name <- baby.names # add the baby names to the data frame
      
      YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
      
    } # End loop over mates
    
  } # End loop over mothers
  
  YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
  
  # Assign each YOY to the cell they are into
  YOY.df.sf <- YOY.df.sf %>% 
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  
  # Rename columns for convenience
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'  
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
  
  
  loopy.pop.sf <- NULL
  loopy.pop.sf <- rbind(init.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
  
  # At end of year 0
  Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  
  # MORTALITY (two steps)
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
  loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                  ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                         ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
  
  # Density-dependent yoy survival
  next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
  diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
  ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
  yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
  loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
  
  # Lethal sampling of individuals. 
  
  # Size of this year's adult & juvenile population (before survival draw).
  number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
  
  # Add the northing (sampling probability) to population dataframe
  loopy.pop.sf$sampling.prob <- sampling.prob[as.numeric(loopy.pop.sf$ID.small),"sampling.prob"]
  
  sample.pop <- loopy.pop.sf %>%
    subset (survival == "M") %>% # Only harvest individuals
    subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
    slice_sample(n = round((harvest.size*number.adults)/100), weight_by = sampling.prob, replace = FALSE) # Harvest n% of this year's total adult population
  
  # Remove northing (I could also just leave the northing to the entire dataset every time..)
  loopy.pop.sf <- subset (loopy.pop.sf, select = - sampling.prob)
  
  # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
  loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                   "H", # Then write him as harvested
                                   loopy.pop.sf$survival) # Otherwise keep S or M as it was before
  
  # Dispersal and grid assignment 
  # Since all my initial pop has dispersed, and none of the YOY for year 0 disperse, there is no dispersal event for the first year
  
  # At end of year 0
  print(paste("year 0 ", "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
  
  # Note: not all males or females are fathers or mothers
  # Mothers and fathers may have contributed to reproduction this year and die afterwards (since survival is estimated at a later stage)
  
  # Visualize the YOY
  # Don't worry if it looks like there are more males than females. The sibling dots overlap, and males are "on top" of females, but sex ratio is okay.
  # ggplot() +
  #   geom_sf(data = sf_box, fill = 'white', lwd = 0.05) +
  #   geom_sf(data = sf_grid, fill = 'transparent', lwd = 0.3) +
  #   geom_sf(data = YOY.df.sf, mapping = aes (color = sex), size = 0.7) #+
  #coord_sf(datum = st_crs(32616)) # To change axis to UTM (meters)
  
  ###---Loop for all other breeding years---###
  
  pop.size <- data.frame()
  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year
  sample.list <- list() # Make list to store dataframe of sampled individual for each year
  
  for(v in 1:(num.years)){ #loop through all of the years in the simulation
    
    v.pop.sf <- loopy.pop.sf[loopy.pop.sf$survival =="S",] # Bring in the data from the previous iteration, but only include those that survive
    v.pop.sf <- subset (v.pop.sf, select = - survival) # Remove the survival column
    v.pop.sf$age <- v.pop.sf$age+1 # Increase each individuals age by one for the new year
    
    mothers <- which(v.pop.sf$sex=='F' & v.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
    
    YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
    
    # Loop through all the available mothers
    for(j in 1:length(mothers)){
      
      num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
      inter.df <- data.frame() # Create a data frame for the offspring born from each mother
      
      # Record the mother
      mother <- as_Spatial(v.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
      
      # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
      fathers <- which(v.pop.sf$sex=='M' & v.pop.sf$age>=repro.age & v.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
      
      if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
      
      # Need to make a combination of number of available fathers and number of mates this year for this female
      ifelse(length(fathers) == 1,  # If only one father available...
             {num.mates.x <- 1; fathers.ID <- v.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
             ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
      
      num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
      num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
      
      if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
      
      # If there are two fathers, need to distribute the cubs among them
      ifelse (num.offspring %% 2 == 0, # If even number of cubs....
              num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
              num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
      
      for (h in 1:num.mates.x){
        
        ifelse (num.mates.x == 1, 
                num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
                ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                       num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                       num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
        
        if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
        
        indv.name <- NA # Create a place holder for the random name given to each offspring
        age <- 0 # Assign a 0 age to each offspring born this year
        birth.year <- v # Note that these pups were born in the 0 year of the simulation
        sex <- NA # Create a place holder for the sexes
        long <- NA
        lat <- NA
        dispersal <- NA
        
        
        inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                     "mother" = mother[[1]], # record the mother's name
                                     "father" = fathers.ID[[1]][[h]], # record the father's name,
                                     sex, long, lat, dispersal) 
        
        inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
        inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
        inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
        inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
        inter.df$dispersal <- 0 # Newborns have not dispersed yet
        
        baby.names <- vector() # Create a blank vector for random baby names
        
        for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
          name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
          baby.names <- c(baby.names,name1)
        } # End loop over name
        
        inter.df$indv.name <- baby.names # add the baby names to the data frame
        
        YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
        
      } # End loop over mates
      
    } # End loop over mothers
    
    YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
    
    # Assign each YOY to the cell they are into
    YOY.df.sf <- YOY.df.sf %>% 
      st_join(sf_grid_small, join = st_intersects)  %>% # Spatial join of each point to the grid ID. Adds a new column.
      st_join(sf_grid_large, join = st_intersects)
    
    
    # Rename columns for convenience
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
    
    loopy.pop.sf <- NULL
    loopy.pop.sf <- rbind(v.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
    
    # At end of year 0
    Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    
    # DISPERSAL
    
    # Let's fix dispersal at 1yo. 
    dispersers.pop <- data.frame() # Create a place holder for this year's dispersers
    
    dispersers.pop <- loopy.pop.sf %>%
      subset (dispersal == 0) %>% # Only one dispersal events, so remove those that have already dispersed (only used if eg only 50% of 1yo dispersed, and the rest at 2yo)
      subset (age == 1) %>% # All 1 year olds disperse
      st_drop_geometry() %>% # Remove geometry for now
      subset (select = - ID.large) # Remove large cell ID (it might change after dispersal)
    
    for (d in 1:nrow(dispersers.pop)) {
      
      dest.vec <- transition.distance.mu[as.numeric(dispersers.pop[d, "ID.small"]),] # Get the cell # of the original location ("from)
      
      dispersers.pop[d, "ID.small"] <- sample(1:nrow(sf_grid_small), size = 1, prob = dest.vec) # Select final destination ("to") based on distance x mu (probability of same point = 0)
      
      # Retrieve geometry from small grid ID (based on centroids), transformed into long and lat columns
      dispersers.pop[d, "long"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "X"]
      dispersers.pop[d, "lat"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "Y"]
    }
    
    dispersers.pop.sf <- st_as_sf(dispersers.pop, coords = c("long", "lat"), crs = 4326) %>%  # Create the spatial feature for dispersers
      st_join(sf_grid_large, join = st_intersects) # Spatial join of each point to the large grid ID. Adds a new column.
    
    # Rename column for convenience
    names(dispersers.pop.sf)[names(dispersers.pop.sf) == 'ID'] <- 'ID.large'
    
    # Assign dispersal as being done
    dispersers.pop.sf$dispersal <- 1
    
    # DIAGNOSTIC: let's connect the "start" and "end" of dispersal. Gives an idea of how far individuals move around.
    #
    # dispersers.pop.origin <- loopy.pop.sf %>% # This one is just for diagnostic purposes
    #   subset (dispersal == 0) %>%
    #   subset (age == 1)
    # 
    # start.long <- st_coordinates(dispersers.pop.origin)[,1]
    # start.lat <- st_coordinates(dispersers.pop.origin)[,2]
    # 
    # end.long <- st_coordinates(dispersers.pop.sf)[,1]
    # end.lat <- st_coordinates(dispersers.pop.sf)[,2]
    # 
    # # plot(c(start.long, end.long), c(start.lat, end.lat), col = rep(c("red", "blue"), each = length(start.long)),
    # #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # # arrows(start.long, start.lat, end.long, end.lat)
    # 
    # # If it's too messy, you can subset the dataframes
    # 
    # start.long.sub <- start.long[1:(length(start.long)/5)]
    # start.lat.sub <- start.lat[1:(length(start.lat)/5)]
    # 
    # end.long.sub <- end.long[1:(length(end.long)/5)]
    # end.lat.sub <- end.lat[1:(length(end.lat)/5)]
    # 
    # plot(c(start.long.sub, end.long.sub), c(start.lat.sub, end.lat.sub), col = rep(c("red", "blue"), each = length(start.long.sub)),
    #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # arrows(start.long.sub, start.lat.sub, end.long.sub, end.lat.sub)
    
    # Remove individuals that are in dispersers.pop from the loopy.pop, and add them again with the new dataset
    loopy.pop.sf <- loopy.pop.sf[!loopy.pop.sf$indv.name %in% dispersers.pop.sf$indv.name,] # Remove
    loopy.pop.sf <- rbind.data.frame(loopy.pop.sf, dispersers.pop.sf) # Add
    
    # MORTALITY (two steps)
    
    #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
    loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                    ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                           ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
    
    # Density-dependent yoy survival
    next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
    diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
    ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
    yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
    loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
    
    # Lethal sampling of individuals. 
    
    # Size of this year's adult & juvenile population (before survival draw).
    number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
    
    # Add the northing (sampling probability) to population dataframe
    loopy.pop.sf$sampling.prob <- sampling.prob[as.numeric(loopy.pop.sf$ID.small),"sampling.prob"]
    
    sample.pop <- loopy.pop.sf %>%
      subset (survival == "M") %>% # Only harvest individuals
      subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
      slice_sample(n = round((harvest.size*number.adults)/100), weight_by = sampling.prob, replace = FALSE) # Harvest n% of this year's total adult population
    
    # Remove northing (I could also just leave the northing to the entire dataset every time..)
    loopy.pop.sf <- subset (loopy.pop.sf, select = - sampling.prob)
    
    # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
    loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                     "H", # Then write him as harvested
                                     loopy.pop.sf$survival) # Otherwise keep S or M as it was before
    
    
    # Return samples to a data frame, with geometry transformed into long and lat columns
    sample.pop <- sample.pop %>%
      mutate(long = st_coordinates(.)[,"X"], #Re-extract long and lat coordinates in separate columns
             lat = st_coordinates(.)[,"Y"]) %>%
      st_set_geometry(NULL) # Remove geometry
    
    sample.list[[v]] <- sample.pop
    loopy.list[[v]] <- loopy.pop.sf # Save the current year's population data as a list element, where the index corresponds to the year
    
    # Print the simulation year and the population size in the R console so they can be observed
    print(paste("year", v, "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df.sf), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
    
    ### Count survivors for each year v
    pop.size.vec <- cbind.data.frame(year=v,
                                     population_size_preharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]), # Minus young of the year
                                     population_size_postharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1 & loopy.pop.sf$survival=="S",]), # Minus young of the year
                                     Nm_preharvest = Nm, # Includes males that may have reproduced this year but were subsequently killed (harvest or natural)
                                     Nf_preharvest = Nf) # Includes females that may have reproduced this year but were subsequently killed (harvest or natural)
    
    pop.size <- rbind(pop.size, pop.size.vec)
    
  } # end loop over sim years
  
  # Label the list elements with the year
  names(loopy.list) <- paste0("year.end.pop.", seq(1:num.years))
  
  return(invisible(list(loopy.list, sample.list, pop.size)))
  
}

IPP.restricted <- function(){
  
  # Most basic population every (aka "baseline scenario", with equal and constant reproductive output)
  # Simple age at maturity, maximum age etc...
  
  ###---Set up your initial population---###
  
  init.pop <- data.frame() # Create a blank data frame which will become the initial population
  
  age <- NULL # Create an empty vector for age
  for(y in 1:length(Nages)) age <- c(age, rep(y, Nages[y])) # Elongate age into population size and fill each row according to ratios in Nages
  
  # Loop that creates the below data for each individual in a population the size of "init.pop.size"
  for(i in 1:init.pop.size){
    
    indv.name <- paste(sample(letters, size = 20, replace = T), collapse="") # generate a random name that is 20 letters long
    mother <- "xxxxx" # The individuals in the initial population do not have known mothers
    father <- "xxxxx" # The individuals in the initial population do not have known fathers
    birth.year <- -1 # this is a place holder for individuals born within the simulation
    sex <- sample(c('F','M'), size = 1, prob = c(init.prop.female, 1-init.prop.female)) # randomly draw sex based on the proportions set in the parameter section
    dispersal <- 1 # Let's assume they all dispersed, no matter the age (only 1 dispersal event in life)
    init.vec <- cbind.data.frame(indv.name, birth.year, age[i], mother, father, sex, dispersal) # Create a row with the attributes for each individual - using cbind.data.frame allows them to keep some numeric and others characters
    init.pop <- rbind(init.pop, init.vec) # place the individual into the init.pop data frame
  }
  
  names(init.pop)[3] <- "age" # Change the name to remove the [i] from the loop
  
  # Have a look at your population
  
  head(init.pop)
  summary(as.factor(init.pop$age)) # Should match Nages
  summary(as.factor(init.pop$sex)) # Should be 50/50
  
  # Create a box (100x100km)
  box <- matrix(c(0, 0, # Create the bounding points
                  0, 1,
                  1, 1,
                  1, 0,
                  0, 0),
                ncol=2, byrow=TRUE)
  
  sf_box <- st_polygon(list(box)) # Join the points to create a polygon
  sf_box <- st_sfc(sf_box, crs = 4326) # Transform polygon into a spatial feature
  
  # Then create two grids on box
  # One larger grid for reproduction
  sf_grid_large  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.1, 0.1)) %>% # Make it 10x10
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326) # Reverts the dataframe to a sf object (now including ID)
  
  # One smaller grid for dispersal and heterogenous point process (/density)
  sf_grid_small  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.05, 0.05)) %>% # Make it 50x50
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326)
  
  # Get centroid for smaller grid
  sf_grid_small_centroid <- st_centroid(sf_grid_small)
  
  # Create random points that will populate the small grid (centroid of each cell)
  points_init <- st_sample(sf_grid_small_centroid, size = nrow(init.pop), replace = TRUE)
  
  # Assign those points to the init.pop dataframe
  init.pop$long <- st_coordinates(points_init)[,1]
  init.pop$lat <- st_coordinates(points_init)[,2]
  
  # Assign each ind to the cell they are into
  init.pop.sf <- st_as_sf (init.pop, coords = c("long", "lat"), crs = 4326) %>% # Transform your dataframe as a sf
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  # Rename columns for convenience
  names(init.pop.sf)[names(init.pop.sf) == 'ID.x'] <- 'ID.small'
  names(init.pop.sf)[names(init.pop.sf) == 'ID.y'] <- 'ID.large'
  
  # Visualize what you just did
  
  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_large, lwd = 0.3) +
    geom_sf(data = sf_grid_small, lwd = 0.3) +
    geom_sf(data = init.pop.sf, mapping = aes (color = sex), size = 0.7)
  
  # Create matrix and fill with distances
  distance.matrix <- st_distance(x = sf_grid_small_centroid, y = sf_grid_small_centroid) # units = m
  distance.matrix[1:10, 1:10] # ~5.5km between adjacent grids
  
  # Transition matrix based on distances only (follows dgamma with shape = 4 and rate = 1.5e-4)
  # Note: might need to do one matrix for males and one matrix for females if dispersal distances are different
  transition.distance = matrix(dgamma(as.numeric(distance.matrix), 2, 1.5e-4), 400, 400)
  
  # Divide by a normalizing constant so that sum(transition.distance) == 1
  transition.distance <- transition.distance / sum(transition.distance)
  
  # Calculate SUM probability of ind moving to each cell based on distance only
  transition.distance.vector <- colSums(transition.distance) 
  sf_grid_small_centroid$TDV <- transition.distance.vector # Using the plot below, you can see that there will naturally be more points in the center of the box
  
  ggplot() +
    geom_sf(data = sf_box,lwd = 0.05) +
    geom_sf(data = sf_grid_small_centroid, mapping = aes (color = TDV))
  
  # Create a density that is the *exact* opposite of the natural density created by distances
  mu_opposite <-  mean(transition.distance.vector) / transition.distance.vector
  
  # Create the intensity parameter mu
  mu <- c (rep(seq(from = 10, to = 1, by = -1), each = 20),
           rep(seq(from = 1, to = 10, by = 1), each = 20))
  
  # Empty matrix
  transition.distance.mu <- matrix(0, 400, 400) 
  
  for (i in 1: nrow(sf_grid_small)) { # For each row ("from")
    for (j in 1:nrow(sf_grid_small)) { # For each column ("to")
      transition.distance.mu[i,j] <- transition.distance[i,j] * mu[j] * mu_opposite[j] # Multiply distance with mu and the mu_opposite
    }
  }
  
  mu <- colSums(transition.distance.mu) # Now you can see that cells in the same row have equal probability of being "moved to", whilst having the IPP
  
  # Create a dataframe for sampling probability // Northing will be sampling probability
  sampling.prob <- data.frame (ID.small = c(1:400), 
                               sampling.prob = c(rep(0, times = 200),
                                                 rep(seq(from = 1, to = 10, by = 1), each = 20)))
  
  ###---Breeding---###
  
  ###---Year 0 breeding---###
  
  mothers <- which(init.pop.sf$sex=='F' & init.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
  
  YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
  
  # Loop through all the available mothers
  for(j in 1:length(mothers)){
    
    num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
    inter.df <- data.frame() # Create a data frame for the offspring born from each mother
    
    # Record the mother
    mother <- as_Spatial(init.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
    
    # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
    fathers <- which(init.pop.sf$sex=='M' & init.pop.sf$age>=repro.age & init.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
    
    if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
    
    # Need to make a combination of number of available fathers and number of mates this year for this female
    ifelse(length(fathers) == 1,  # If only one father available...
           {num.mates.x <- 1; fathers.ID <- init.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
           ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                   fathers.ID <- init.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
    
    num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
    num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
    
    if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
    
    # If there are two fathers, need to distribute the cubs among them
    ifelse (num.offspring %% 2 == 0, # If even number of cubs....
            num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
            num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
    
    for (h in 1:num.mates.x){
      
      ifelse (num.mates.x == 1, 
              num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
              ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                     num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                     num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
      
      if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
      
      indv.name <- NA # Create a place holder for the random name given to each offspring
      age <- 0 # Assign a 0 age to each offspring born this year
      birth.year <- 0 # Note that these pups were born in the 0 year of the simulation
      sex <- NA # Create a place holder for the sexes
      long <- NA
      lat <- NA
      dispersal <- NA
      
      
      inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                   "mother" = mother[[1]], # record the mother's name
                                   "father" = fathers.ID[[1]][[h]], # record the father's name,
                                   sex, long, lat, dispersal) 
      
      inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
      inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
      inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
      inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
      inter.df$dispersal <- 0 # Newborns have not dispersed yet
      
      baby.names <- vector() # Create a blank vector for random baby names
      
      for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
        name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
        baby.names <- c(baby.names,name1)
      } # End loop over name
      
      inter.df$indv.name <- baby.names # add the baby names to the data frame
      
      YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
      
    } # End loop over mates
    
  } # End loop over mothers
  
  YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
  
  # Assign each YOY to the cell they are into
  YOY.df.sf <- YOY.df.sf %>% 
    st_join(sf_grid_small, join = st_intersects) %>% # Spatial join of each point to the grid ID. Adds a new column.
    st_join(sf_grid_large, join = st_intersects) 
  
  
  # Rename columns for convenience
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'  
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
  
  
  loopy.pop.sf <- NULL
  loopy.pop.sf <- rbind(init.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
  
  # At end of year 0
  Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
  
  # MORTALITY (two steps)
  
  #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
  loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                  ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                         ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
  
  # Density-dependent yoy survival
  next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
  diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
  ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
  yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
  loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
  
  # Lethal sampling of individuals. 
  
  # Size of this year's adult & juvenile population (before survival draw).
  number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
  
  # Add the northing (sampling probability) to population dataframe
  loopy.pop.sf$sampling.prob <- sampling.prob[as.numeric(loopy.pop.sf$ID.small),"sampling.prob"]
  
  sample.pop <- loopy.pop.sf %>%
    subset (survival == "M") %>% # Only harvest individuals
    subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
    slice_sample(n = round((harvest.size*number.adults)/100), weight_by = sampling.prob, replace = FALSE) # Harvest n% of this year's total adult population
  
  # Remove northing (I could also just leave the northing to the entire dataset every time..)
  loopy.pop.sf <- subset (loopy.pop.sf, select = - sampling.prob)
  
  # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
  loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                   "H", # Then write him as harvested
                                   loopy.pop.sf$survival) # Otherwise keep S or M as it was before
  
  # Dispersal and grid assignment 
  # Since all my initial pop has dispersed, and none of the YOY for year 0 disperse, there is no dispersal event for the first year
  
  # At end of year 0
  print(paste("year 0 ", "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
  
  # Note: not all males or females are fathers or mothers
  # Mothers and fathers may have contributed to reproduction this year and die afterwards (since survival is estimated at a later stage)
  
  # Visualize the YOY
  # Don't worry if it looks like there are more males than females. The sibling dots overlap, and males are "on top" of females, but sex ratio is okay.
  # ggplot() +
  #   geom_sf(data = sf_box, fill = 'white', lwd = 0.05) +
  #   geom_sf(data = sf_grid, fill = 'transparent', lwd = 0.3) +
  #   geom_sf(data = YOY.df.sf, mapping = aes (color = sex), size = 0.7) #+
  #coord_sf(datum = st_crs(32616)) # To change axis to UTM (meters)
  
  ###---Loop for all other breeding years---###
  
  pop.size <- data.frame()
  loopy.list <- list() # Make list to store dataframe of population for each year, where each element corresponds to the year e.g. loopy.list[[1]] is the population from the first year
  sample.list <- list() # Make list to store dataframe of sampled individual for each year
  
  for(v in 1:(num.years)){ #loop through all of the years in the simulation
    
    v.pop.sf <- loopy.pop.sf[loopy.pop.sf$survival =="S",] # Bring in the data from the previous iteration, but only include those that survive
    v.pop.sf <- subset (v.pop.sf, select = - survival) # Remove the survival column
    v.pop.sf$age <- v.pop.sf$age+1 # Increase each individuals age by one for the new year
    
    mothers <- which(v.pop.sf$sex=='F' & v.pop.sf$age>=repro.age) # Determine which females are available to breed in this year
    
    YOY.df <- data.frame() # Create an empty data frame to populate with the YOY born in this year
    
    # Loop through all the available mothers
    for(j in 1:length(mothers)){
      
      num.mates.x <- sample(num.mates, size = 1) # Determine the number of mates for a given mother this year
      inter.df <- data.frame() # Create a data frame for the offspring born from each mother
      
      # Record the mother
      mother <- as_Spatial(v.pop.sf[mothers[j],]) # Note: the geometry is "attached" to the row if working with sf
      
      # Determine which fathers are available to breed in this year (for each specifically female due to spatial stuff)
      fathers <- which(v.pop.sf$sex=='M' & v.pop.sf$age>=repro.age & v.pop.sf$ID.large == mother$ID.large) # All mature males that are within the same grid as the female
      
      if(length(fathers)==0){next} # if there are no available father in the grid (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
      
      # Need to make a combination of number of available fathers and number of mates this year for this female
      ifelse(length(fathers) == 1,  # If only one father available...
             {num.mates.x <- 1; fathers.ID <- v.pop.sf[fathers, 1]}, # ...he's the only fatherand save his name (multiple actions on the TRUE statement)
             ifelse( num.mates.x == 1, # If more than one father available, but only one mate for this female...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 1, replace = FALSE), 1], #... randomly select one male...
                     fathers.ID <- v.pop.sf[sample(fathers, size = 2, replace = FALSE), 1])) #... otherwise select two.
      
      num.offspring <- round(rnorm(1, f, 0.5)) # Generate litter size
      num.offspring[which(num.offspring < 0)] = 0 # Ensure that there are no negative values for num.offspring (otherwise it breaks the loop)
      
      if(num.offspring==0){next} # if there were no offspring from this mother, skips to the next mother (if you don't include this you will get an error in the loop)
      
      # If there are two fathers, need to distribute the cubs among them
      ifelse (num.offspring %% 2 == 0, # If even number of cubs....
              num.offspring.even <- (num.offspring / 2), # ... divide evenly between the two parents.
              num.offspring.odd <- c (ceiling(num.offspring / 2) , floor(num.offspring / 2)))  # If odd litter size, distribute unevenly. Round up with ceiling, and round down with floor. Both functions round to the nearest integer but in a different direction.
      
      for (h in 1:num.mates.x){
        
        ifelse (num.mates.x == 1, 
                num.offspring.mate <- num.offspring,  # If only 1 mate, all the offspring belong to him
                ifelse(num.offspring %% 2 == 0, # If even number of cubs....
                       num.offspring.mate <- num.offspring.even, # ... even number of cubs for each mate ...
                       num.offspring.mate <- num.offspring.odd[h])) # ... else distribute based on odd or even number of cubs
        
        if(num.offspring.mate==0){next} # if there were no offspring from this pair, skips to the next mother (if you don't include this you will get an error in the loop)
        
        indv.name <- NA # Create a place holder for the random name given to each offspring
        age <- 0 # Assign a 0 age to each offspring born this year
        birth.year <- v # Note that these pups were born in the 0 year of the simulation
        sex <- NA # Create a place holder for the sexes
        long <- NA
        lat <- NA
        dispersal <- NA
        
        
        inter.df <- cbind.data.frame(indv.name, birth.year, age, # Create a row with the attributes of each YOY (at this point only one from each mating pair)
                                     "mother" = mother[[1]], # record the mother's name
                                     "father" = fathers.ID[[1]][[h]], # record the father's name,
                                     sex, long, lat, dispersal) 
        
        inter.df <- inter.df[rep(seq_len(nrow(inter.df)), num.offspring.mate), ] # Create the random number of offspring from each pairing that we set above
        inter.df$sex <- sample(c('F','M'), size = nrow(inter.df), prob = birth.sex.ratio, replace = T) # Assign biological sex to each new born based on the sex ration set in the parameter section
        inter.df$long <- st_coordinates(st_as_sf(mother))[,1] # Give the mother's coordinates to the YOY
        inter.df$lat <- st_coordinates(st_as_sf(mother))[,2]
        inter.df$dispersal <- 0 # Newborns have not dispersed yet
        
        baby.names <- vector() # Create a blank vector for random baby names
        
        for(w in 1:nrow(inter.df)){ # create a random name for each newborn from this mating pair
          name1 <- paste(sample(letters, size = 20, replace = T), collapse="")
          baby.names <- c(baby.names,name1)
        } # End loop over name
        
        inter.df$indv.name <- baby.names # add the baby names to the data frame
        
        YOY.df <- rbind(YOY.df,inter.df) # Add all of the offspring from this mother/mate pairing to a data frame of all the YOY for the year
        
      } # End loop over mates
      
    } # End loop over mothers
    
    YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 4326) # Transform your dataframe as a sf to bind with existing population
    
    # Assign each YOY to the cell they are into
    YOY.df.sf <- YOY.df.sf %>% 
      st_join(sf_grid_small, join = st_intersects)  %>% # Spatial join of each point to the grid ID. Adds a new column.
      st_join(sf_grid_large, join = st_intersects)
    
    
    # Rename columns for convenience
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.x'] <- 'ID.small'
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID.y'] <- 'ID.large'
    
    loopy.pop.sf <- NULL
    loopy.pop.sf <- rbind(v.pop.sf, YOY.df.sf) #Combine the YOY data with the other individuals present this year
    
    # At end of year 0
    Nf <- sum(loopy.pop.sf$sex=='F' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    Nm <- sum(loopy.pop.sf$sex=='M' & loopy.pop.sf$age >= repro.age, na.rm=TRUE)
    
    # DISPERSAL
    
    # Let's fix dispersal at 1yo. 
    dispersers.pop <- data.frame() # Create a place holder for this year's dispersers
    
    dispersers.pop <- loopy.pop.sf %>%
      subset (dispersal == 0) %>% # Only one dispersal events, so remove those that have already dispersed (only used if eg only 50% of 1yo dispersed, and the rest at 2yo)
      subset (age == 1) %>% # All 1 year olds disperse
      st_drop_geometry() %>% # Remove geometry for now
      subset (select = - ID.large) # Remove large cell ID (it might change after dispersal)
    
    for (d in 1:nrow(dispersers.pop)) {
      
      dest.vec <- transition.distance.mu[as.numeric(dispersers.pop[d, "ID.small"]),] # Get the cell # of the original location ("from)
      
      dispersers.pop[d, "ID.small"] <- sample(1:nrow(sf_grid_small), size = 1, prob = dest.vec) # Select final destination ("to") based on distance x mu (probability of same point = 0)
      
      # Retrieve geometry from small grid ID (based on centroids), transformed into long and lat columns
      dispersers.pop[d, "long"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "X"]
      dispersers.pop[d, "lat"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "Y"]
    }
    
    dispersers.pop.sf <- st_as_sf(dispersers.pop, coords = c("long", "lat"), crs = 4326) %>%  # Create the spatial feature for dispersers
      st_join(sf_grid_large, join = st_intersects) # Spatial join of each point to the large grid ID. Adds a new column.
    
    # Rename column for convenience
    names(dispersers.pop.sf)[names(dispersers.pop.sf) == 'ID'] <- 'ID.large'
    
    # Assign dispersal as being done
    dispersers.pop.sf$dispersal <- 1
    
    # DIAGNOSTIC: let's connect the "start" and "end" of dispersal. Gives an idea of how far individuals move around.
    #
    # dispersers.pop.origin <- loopy.pop.sf %>% # This one is just for diagnostic purposes
    #   subset (dispersal == 0) %>%
    #   subset (age == 1)
    # 
    # start.long <- st_coordinates(dispersers.pop.origin)[,1]
    # start.lat <- st_coordinates(dispersers.pop.origin)[,2]
    # 
    # end.long <- st_coordinates(dispersers.pop.sf)[,1]
    # end.lat <- st_coordinates(dispersers.pop.sf)[,2]
    # 
    # # plot(c(start.long, end.long), c(start.lat, end.lat), col = rep(c("red", "blue"), each = length(start.long)),
    # #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # # arrows(start.long, start.lat, end.long, end.lat)
    # 
    # # If it's too messy, you can subset the dataframes
    # 
    # start.long.sub <- start.long[1:(length(start.long)/5)]
    # start.lat.sub <- start.lat[1:(length(start.lat)/5)]
    # 
    # end.long.sub <- end.long[1:(length(end.long)/5)]
    # end.lat.sub <- end.lat[1:(length(end.lat)/5)]
    # 
    # plot(c(start.long.sub, end.long.sub), c(start.lat.sub, end.lat.sub), col = rep(c("red", "blue"), each = length(start.long.sub)),
    #      pch = 20, cex = 2, xlab = "lon", ylab = "lat")
    # arrows(start.long.sub, start.lat.sub, end.long.sub, end.lat.sub)
    
    # Remove individuals that are in dispersers.pop from the loopy.pop, and add them again with the new dataset
    loopy.pop.sf <- loopy.pop.sf[!loopy.pop.sf$indv.name %in% dispersers.pop.sf$indv.name,] # Remove
    loopy.pop.sf <- rbind.data.frame(loopy.pop.sf, dispersers.pop.sf) # Add
    
    # MORTALITY (two steps)
    
    #Assign a column with whether or not an individual this year will survive to next year. The survival is randomly drawn based on the life stage and the probabilities defined in the parameter section
    loopy.pop.sf$survival <- ifelse(loopy.pop.sf$age==0, "S",
                                    ifelse(loopy.pop.sf$age<repro.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age > 0 & loopy.pop.sf$age<repro.age)), prob=c(juvenile.survival, 1-juvenile.survival),replace=T),
                                           ifelse(loopy.pop.sf$age<max.age, sample(c("S","M"), size=length(which(loopy.pop.sf$age >= repro.age & loopy.pop.sf$age<max.age)), prob=c(adult.survival, 1-adult.survival),replace=T), "M"))) 
    
    # Density-dependent yoy survival
    next.pop.size <- as.numeric(count(loopy.pop.sf[loopy.pop.sf$survival == "S",])) # next year's population size if all yoy survived
    diff.pop.size <- next.pop.size - init.pop.size # This is how many yoy must die
    ratio <- diff.pop.size / as.numeric(nrow(YOY.df)) # proportion of yoy that must die
    yoy.survival <- rnorm(1, mean = 1- ratio, sd = 0.01) # This year's yoy survival, drawn from Normal distribution. Add just a small amount of stochasticity so that the population is not **always** 3,000?
    loopy.pop.sf$survival[loopy.pop.sf$age==0] <- sample(c("S", "M"), size = length (which(loopy.pop.sf$age == 0)), prob = c(yoy.survival, 1-yoy.survival), replace = T) # Assign survival to YOYs
    
    # Lethal sampling of individuals. 
    
    # Size of this year's adult & juvenile population (before survival draw).
    number.adults = nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]) 
    
    # Add the northing (sampling probability) to population dataframe
    loopy.pop.sf$sampling.prob <- sampling.prob[as.numeric(loopy.pop.sf$ID.small),"sampling.prob"]
    
    sample.pop <- loopy.pop.sf %>%
      subset (survival == "M") %>% # Only harvest individuals
      subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
      slice_sample(n = round((harvest.size*number.adults)/100), weight_by = sampling.prob, replace = FALSE) # Harvest n% of this year's total adult population
    
    # Remove northing (I could also just leave the northing to the entire dataset every time..)
    loopy.pop.sf <- subset (loopy.pop.sf, select = - sampling.prob)
    
    # Now, assign the sampled individuals to Survival = H (harvested) in the loopy.pop.sf dataframe
    loopy.pop.sf$survival <- ifelse (loopy.pop.sf$indv.name %in% sample.pop$indv.name, # If ind.name is in samplepop
                                     "H", # Then write him as harvested
                                     loopy.pop.sf$survival) # Otherwise keep S or M as it was before
    
    
    # Return samples to a data frame, with geometry transformed into long and lat columns
    sample.pop <- sample.pop %>%
      mutate(long = st_coordinates(.)[,"X"], #Re-extract long and lat coordinates in separate columns
             lat = st_coordinates(.)[,"Y"]) %>%
      st_set_geometry(NULL) # Remove geometry
    
    sample.list[[v]] <- sample.pop
    loopy.list[[v]] <- loopy.pop.sf # Save the current year's population data as a list element, where the index corresponds to the year
    
    # Print the simulation year and the population size in the R console so they can be observed
    print(paste("year", v, "Nf=", Nf, "Nm=", Nm, "N_newborns=", nrow(YOY.df.sf), "N_deaths=", sum(loopy.pop.sf$survival=="M", loopy.pop.sf$survival=="H"), "N_survivors= ", nrow(loopy.pop.sf[loopy.pop.sf$survival=="S",]) , sep=" "))
    
    ### Count survivors for each year v
    pop.size.vec <- cbind.data.frame(year=v,
                                     population_size_preharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1,]), # Minus young of the year
                                     population_size_postharvest=nrow(loopy.pop.sf[loopy.pop.sf$age >= 1 & loopy.pop.sf$survival=="S",]), # Minus young of the year
                                     Nm_preharvest = Nm, # Includes males that may have reproduced this year but were subsequently killed (harvest or natural)
                                     Nf_preharvest = Nf) # Includes females that may have reproduced this year but were subsequently killed (harvest or natural)
    
    pop.size <- rbind(pop.size, pop.size.vec)
    
  } # end loop over sim years
  
  # Label the list elements with the year
  names(loopy.list) <- paste0("year.end.pop.", seq(1:num.years))
  
  return(invisible(list(loopy.list, sample.list, pop.size)))
  
}
