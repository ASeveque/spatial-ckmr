####---------- Load packages and source functions ----------####

library(tidyverse) # Data management
library(sf) # Spatial stuff
library(terra) # Spatial stuff (for rasters)
library(geosphere) # Spatial stuff (for simulations only)
library(ggplot2) # For plots
library(viridis)

####---------- Demographic parameters ----------####

# Set a seed
set.seed(123456)

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
init.adult.pop.size <- 100 # Initial number of mothers

# Create an age distribution based on average fecundity and survival
f <- 3 # Average annual adult fecundity (= number of offspring)
props <- rep (NA, max.age+1) # Placeholder
props[1] <- f # YOY (those just being born)
props[2] <- f * yoy.survival # YOY survival (from 0 to 1 years)
for (y in 3:(repro.age+1)) props[y] <- props[y-1] * juvenile.survival # Juvenile survival (from 1 to adult years)
for (y in (repro.age+2):(max.age+1)) props[y] <- props[y-1] * adult.survival # Adult survival

Nages <- round(props[-1] * init.adult.pop.size) # Transforms age "ratios" to number of individuals in each age + remove YOY
init.pop.size <- sum(Nages) # Total initial population size

# Sampling variables
harvest.size = 10 # Harvest n% of the previous year's population
harvest.age = 1 # Minimum age of harvest
sample_first_y = num.years - 2 # First year from simulation where we start genetic sampling
sample_subset <- seq(from=sample_first_y, to=num.years, by=1) # Years of simulation where we will keep samples for CKMR. "by=n" is the sampling frequency.
estimation_year = 28

# Load Michigan mainland
UP <- st_read("./Data/UP_Michigan_mainland.shp") %>% # Upper Peninsula borders
  st_transform(32616)

# Load Michigan counties
counties <- st_read("./Data/Counties_(v17a).shp") %>% 
  st_transform(32616) %>%
  st_intersection(UP) 

# Find 4 biggest counties
polygons_df <- counties %>%
  mutate(area = as.numeric(st_area(.))) %>%
  st_drop_geometry() %>%
  bind_cols(geometry = st_geometry(counties))

largest_four <- polygons_df %>%
  arrange(desc(area)) %>%
  slice_head(n = 4)

largest_four_sf <- st_as_sf(largest_four)

# Look
# ggplot() +
#   geom_sf(data = UP) +
#   geom_sf(data = counties) +
#   geom_sf(data = largest_four_sf, fill = "red")


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
  
  # One smaller grid for dispersal and heterogenous point process (/density)
  sf_grid_small  <- st_make_grid(UP, # Grid will follow the box's bounding box and crs (if any)
                                 square = FALSE,
                                 cellsize = c(10000, 10000))%>% 
    st_intersection(UP) %>%
    st_sf(crs = 32616)
  
  # Calculate the area of each sf_grid_small cell
  intersected_areas <- st_area(sf_grid_small)
  
  # Calculate the area of a full sf_grid_small cell
  full_cell_area <- 10000 * 10000
  
  # Calculate the percentage of each cell that's within UP
  sf_grid_small$percentage_within <- as.numeric(intersected_areas / full_cell_area)
  
  # Keep only cells that are more than 50% within UP (this will remove all the annoying bits on the edge)
  sf_grid_small <- sf_grid_small[sf_grid_small$percentage_within > 0.5,]
  
  # Remove percentage_within
  sf_grid_small <- subset(sf_grid_small, select = -c(percentage_within))
 
  # Add ID column
  sf_grid_small$ID <- 1:nrow(sf_grid_small)
  
  # Get centroid for smaller grid
  sf_grid_small_centroid <- st_centroid(sf_grid_small)
  
  # Create random points that will populate the small grid (centroid of each cell)
  points_init <- st_sample(sf_grid_small_centroid, size = nrow(init.pop), replace = TRUE)
  
  # Assign those points to the init.pop dataframe
  init.pop$long <- st_coordinates(points_init)[,1]
  init.pop$lat <- st_coordinates(points_init)[,2]
  
  # Assign each ind to the cell they are into
  init.pop.sf <- st_as_sf (init.pop, coords = c("long", "lat"), crs = 32616) %>% # Transform your dataframe as a sf
    st_join(sf_grid_small, join = st_intersects) # Spatial join of each point to the grid ID. Adds a new column.

  # Rename columns for convenience
  names(init.pop.sf)[names(init.pop.sf) == 'ID'] <- 'ID.small'

  # Visualize what you just did
  
  # ggplot() +
  #   geom_sf(data = UP,lwd = 0.05) +
  #   geom_sf(data = sf_grid_small, lwd = 0.3) +
  #   geom_sf(data = init.pop.sf, mapping = aes (color = sex), size = 0.7)
  
  # Create matrix and fill with distances
  distance.matrix <- st_distance(x = sf_grid_small_centroid, y = sf_grid_small_centroid) # units = m
  distance.matrix[1:10, 1:10] 
  
  # Transition matrices based on distance only
  transition.distance.female = matrix(dgamma(as.numeric(distance.matrix), 1.5, 0.1e-3), nrow(distance.matrix), nrow(distance.matrix))
  transition.distance.male = matrix(dgamma(as.numeric(distance.matrix), 5, 1e-4), nrow(distance.matrix), nrow(distance.matrix))
  
  # Divide by a normalizing constant so that sum(transition.distance) == 1
  transition.distance.female <- transition.distance.female / sum(transition.distance.female)
  transition.distance.male <- transition.distance.male / sum(transition.distance.male)
  
  # Create the intensity parameter mu
  
  # Let's make a population that is 2x as populated to the west than 
  # it i to the east 
  pop <- sf_grid_small_centroid %>%
    st_transform(crs = 4326)
  pop$long <- st_coordinates(pop)[,1]
  mu <- ifelse(pop$long <= median(pop$long), 2, 1)
  
  # Empty matrix
  transition.distance.female.mu <- matrix(0, nrow(distance.matrix), nrow(distance.matrix)) 
  transition.distance.male.mu <- matrix(0, nrow(distance.matrix), nrow(distance.matrix)) 
  
  for (i in 1: nrow(sf_grid_small)) { # For each row ("from")
    for (j in 1:nrow(sf_grid_small)) { # For each column ("to")
      transition.distance.female.mu[i,j] <- transition.distance.female[i,j] * mu[j]  # Multiply distance with mu
    }
  }
  
  for (i in 1: nrow(sf_grid_small)) { # For each row ("from")
    for (j in 1:nrow(sf_grid_small)) { # For each column ("to")
      transition.distance.male.mu[i,j] <- transition.distance.male[i,j] * mu[j]  # Multiply distance with mu
    }
  }
  
  # Calculate SUM probability of ind moving to each cell based on distance only
  transition.distance.female.vector <- colSums(transition.distance.female.mu) 
  sf_grid_small_centroid$TDV.female <- transition.distance.female.vector 
  transition.distance.male.vector <- colSums(transition.distance.male.mu) 
  sf_grid_small_centroid$TDV.male <- transition.distance.male.vector 
  
  # Check
  # ggplot() +
  #   geom_sf(data = sf_grid_small_centroid, aes(col = TDV.male)) +
  #   scale_color_gradient(low = "green", high = "red")
  
  # Create a dataframe for sampling probability
  # Only cells that are within the largest_four_sf counties are sampled
  intersections <- st_intersects(sf_grid_small_centroid, largest_four_sf)
  sf_grid_small_centroid$sampled <- sapply(intersections, function(x) ifelse(length(x) > 0, 1, 0))
  
  # Separate dataframe
  sampling.prob <- data.frame (ID.small = c(1:nrow(distance.matrix)), 
                               sampling.prob = sf_grid_small_centroid$sampled)
  
  # Check
  # ggplot() +
  #   geom_sf(data = sf_grid_small_centroid, aes(col = sampled))
  
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
    
    # Determine which fathers are available to breed in this year
    fathers <- which(init.pop.sf$sex=='M' & init.pop.sf$age>=repro.age)
    
    if(length(fathers)==0){next} # if there are no available father (somehow), skips to the next mating pair (if you don't include this you will get an error in the loop)
    
    # Need to make a combination of number of available fathers and number of mates this year for this female
    ifelse(length(fathers) == 1,  # If only one father available...
           {num.mates.x <- 1; fathers.ID <- init.pop.sf[fathers, 1]}, # ...he's the only father and save his name (multiple actions on the TRUE statement)
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
  
  YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 32616) # Transform your dataframe as a sf to bind with existing population
  
  # Assign each YOY to the cell they are into
  YOY.df.sf <- YOY.df.sf %>% 
    st_join(sf_grid_small, join = st_intersects) # Spatial join of each point to the grid ID. Adds a new column.

  # Rename columns for convenience
  names(YOY.df.sf)[names(YOY.df.sf) == 'ID'] <- 'ID.small'  

  loopy.pop.sf <- NULL
  loopy.pop.sf <- rbind(init.pop.sf, YOY.df.sf) # Combine the YOY data with the other individuals present this year
  
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
  
  # Add the sampling.prob (sampling probability) to population dataframe
  loopy.pop.sf$sampling.prob <- sampling.prob[as.numeric(loopy.pop.sf$ID.small),"sampling.prob"]
  
  # Normalizing constant 
  loopy.pop.sf$sampling.prob <- loopy.pop.sf$sampling.prob / sum(loopy.pop.sf$sampling.prob)
  
  sample.pop <- loopy.pop.sf %>%
    subset (survival == "M") %>% # Only harvest individuals
    subset (age >= 1) %>% # Only harvest juveniles and adults (no YOY)
    slice_sample(n = round((harvest.size*number.adults)/100), weight_by = sampling.prob, replace = FALSE) # Harvest n% of this year's total adult population
  
  # Remove sampling.prob (I could also just leave the sampling.prob to the entire dataset every time..)
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
      fathers <- which(v.pop.sf$sex=='M' & v.pop.sf$age>=repro.age) 
      
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
    
    YOY.df.sf <- st_as_sf (YOY.df, coords = c("long", "lat"), crs = 32616) # Transform your dataframe as a sf to bind with existing population
    
    # Assign each YOY to the cell they are into
    YOY.df.sf <- YOY.df.sf %>% 
      st_join(sf_grid_small, join = st_intersects) # Spatial join of each point to the grid ID. Adds a new column.

    
    # Rename columns for convenience
    names(YOY.df.sf)[names(YOY.df.sf) == 'ID'] <- 'ID.small'

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
      st_drop_geometry() # Remove geometry for now


    for (d in 1:nrow(dispersers.pop)) {
      
      if(dispersers.pop[d,"sex"] == "F") {
      
      dest.vec <- transition.distance.female.mu[as.numeric(dispersers.pop[d, "ID.small"]),] # Get the cell # of the original location ("from)
      
      dest.vec <- dest.vec / sum(dest.vec) # Normalizing constant
      
      dispersers.pop[d, "ID.small"] <- sample(1:nrow(sf_grid_small), size = 1, prob = dest.vec) # Select final destination ("to") based on distance x mu (probability of same point = 0)
      
      # Retrieve geometry from small grid ID (based on centroids), transformed into long and lat columns
      dispersers.pop[d, "long"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "X"]
      dispersers.pop[d, "lat"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "Y"]
      
      }
      
      if(dispersers.pop[d,"sex"] == "M") {
        
        dest.vec <- transition.distance.male.mu[as.numeric(dispersers.pop[d, "ID.small"]),] # Get the cell # of the original location ("from)
        
        dest.vec <- dest.vec / sum(dest.vec) # Normalizing constant
        
        dispersers.pop[d, "ID.small"] <- sample(1:nrow(sf_grid_small), size = 1, prob = dest.vec) # Select final destination ("to") based on distance x mu (probability of same point = 0)
        
        # Retrieve geometry from small grid ID (based on centroids), transformed into long and lat columns
        dispersers.pop[d, "long"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "X"]
        dispersers.pop[d, "lat"] <- st_coordinates(sf_grid_small_centroid)[as.numeric(dispersers.pop[d, "ID.small"]), "Y"]
        
      }
      
    }
    
    dispersers.pop.sf <- st_as_sf(dispersers.pop, coords = c("long", "lat"), crs = 32616)# Create the spatial feature for dispersers

    # Assign dispersal as being done
    dispersers.pop.sf$dispersal <- 1
    
    # Remove individuals that are in dispersers.pop from the loopy.pop, and add them again with the new dataset
    loopy.pop.sf <- loopy.pop.sf[!loopy.pop.sf$indv.name %in% dispersers.pop.sf$indv.name,] # Remove
    loopy.pop.sf <- rbind.data.frame(loopy.pop.sf, dispersers.pop.sf) # Add
    
    loopy.pop.sf <- loopy.pop.sf[!is.na(loopy.pop.sf$ID.small),]

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
    
    # Normalizing constant 
    loopy.pop.sf$sampling.prob <- loopy.pop.sf$sampling.prob / sum(loopy.pop.sf$sampling.prob)
    
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
  
  # Select your samples that will be used to run the CKMR.
  
  samples <- data.frame() # Placeholder
  
  for(s in sample_subset){ # For each year of sampling
    samples_year <- sample.list[[s]] # Retrieve sampled individuals
    samples_year$sampling_year <- rep(s) # Keep track of what year they were sampled
    samples <- rbind(samples, samples_year) # Bind
  }
  
  # Sort dataframe by birth year, so that individual 1 is always older than individual 2 (because it comes first in the df)
  samples <- samples[order(samples$birth.year),]
  
  # Change the names to fit the functions
  colnames(samples)[1] <- "ID" 
  colnames(samples)[2] <- "birth.year" 
  colnames(samples)[3] <- "sampling.age" 
  colnames(samples)[6] <- "sex" 
  colnames(samples)[11] <- "lon" 
  colnames(samples)[13] <- "sampling.year" 

   # Get coordinates in long lat
  samples.sf <- st_as_sf(samples, coords = c("lon", "lat"), crs = 32616) %>%
    st_transform(crs = 4326)
  
  samples$lon <- st_coordinates(samples.sf)[,1]
  samples$lat <- st_coordinates(samples.sf)[,2]
  
  
  # Change the name of the sampled individuals because right now it is a bit strange
  list_names <- unique(c(samples$ID, samples$mother, samples$father))
  new_names <- paste0("sample_", 1:length(list_names))
  
  # For ID column
  df <- data.frame(ID = list_names, 
                   ID_new = new_names)
    result_df <- samples %>%
    left_join(df, by = "ID") %>%
    mutate(ID = coalesce(ID_new, ID)) %>%
    select(-ID_new)
  
  df <- data.frame(mother = list_names, 
                   mother_new = new_names)
    result_df <- result_df %>%
    left_join(df, by = "mother") %>%
    mutate(mother = coalesce(mother_new, mother)) %>%
    select(-mother_new)
  
   df <- data.frame(father = list_names, 
                    father_new = new_names)
   result_df <- result_df %>%
    left_join(df, by = "father") %>%
    mutate(father = coalesce(father_new, father)) %>%
    select(-father_new)
    
   samples <- result_df  
  
  # ggplot() + 
  #   geom_sf(data = samples.sf)
  
  # Run pairwise comparisons just to get MOPs
  
  # Create dataframe of pairwise comparisons with just individual IDs
  
  pairwise.df <- data.frame(t(combn(samples$ID, m=2))) # Generates all the possible combinations. m = 2 means pair comparisons. t(x) gives the dimensions of the matrix (dataframe).
  colnames(pairwise.df) <- c("Ind_1", "Ind_2") # Rename columns so they can easily be joined
 
   # Create dataframe that will be used to extract the birth years for the younger individual from each pairwise comparison using joins.
  Ind1_birth.years <- samples[,c("ID", "birth.year", "sampling.year", "mother", "father", "sex", "lon", "lat")] # Select relevant columns only
  colnames(Ind1_birth.years) <- c("Ind_1","Ind_1_birth","Ind_1_sampling", "Ind_1_mom", "Ind_1_dad", "Ind_1_sex", "Ind_1_lon", "Ind_1_lat") # Rename columns (make sure they are in the same order as the previous command)
  
  Ind2_birth.years <- samples[,c("ID", "birth.year", "sampling.year", "mother", "father", "sex", "lon", "lat")] # Select relevant columns only
  colnames(Ind2_birth.years) <- c("Ind_2","Ind_2_birth","Ind_2_sampling", "Ind_2_mom", "Ind_2_dad", "Ind_2_sex", "Ind_2_lon", "Ind_2_lat") # Rename columns (make sure they are in the same order as the previous command)
  
  # Combine the two dataframes above to extract birth year and parents for each individual in the pairwise comparison matrix. 
  # This is the main pairwise comparison matrix with all comparisons (not all relevant for now) and individual data.
  
  pairwise.df_all <- merge(x = pairwise.df, y = Ind1_birth.years, 
                           by = "Ind_1", sort = FALSE)
  
  pairwise.df_all <- merge(x = pairwise.df_all, y = Ind2_birth.years, 
                           by = "Ind_2", sort = FALSE)
  
  # Change order of columns (for ease of read and logic)
  pairwise.df_all <- pairwise.df_all[, c("Ind_1", "Ind_1_sex", "Ind_1_birth", "Ind_1_sampling", "Ind_1_mom", "Ind_1_dad", "Ind_1_lon", "Ind_1_lat",
                                         "Ind_2", "Ind_2_sex", "Ind_2_birth", "Ind_2_sampling", "Ind_2_mom", "Ind_2_dad", "Ind_2_lon", "Ind_2_lat")]
  
  pop.pairwise.df_filt <- pairwise.df_all %>%
    filter ((Ind_1_sampling - Ind_1_birth) >= repro.age) %>% # Remove potential parents that were harvested before reaching Age at maturity
    filter ((Ind_2_birth - Ind_1_birth) >= repro.age) %>% # Remove potential parents that did not reach Age at maturity when offspring was born
    filter ((Ind_2_birth - Ind_1_birth) <= max.age) %>% # Remove potential parents that were dead when offspring was born (can happen if more than one sampling occasion)
    filter (Ind_1_sampling >= Ind_2_birth) %>% # Remove potential parents that were sampled before offspring was born (because of lethal sampling). For now, let's assume that adults can reproduce and then die in the same year.
    filter(Ind_1_birth != Ind_2_birth) # Filter intra-cohort comparisons (redundant with line 2 but keep it for clarity's sake)
  
  # Separate mothers
  pop.pairwise.df_filt_female <- pop.pairwise.df_filt[pop.pairwise.df_filt$Ind_1_sex == "F",]
  
  # Assign kinships
  pop.pairwise.df_filt_female$KinPair <- ifelse(pop.pairwise.df_filt_female$Ind_1 == pop.pairwise.df_filt_female$Ind_2_mom & pop.pairwise.df_filt_female$Ind_2_sex == "F", "MomDau", # Mother-daughter
                                                ifelse(pop.pairwise.df_filt_female$Ind_1 == pop.pairwise.df_filt_female$Ind_2_mom & pop.pairwise.df_filt_female$Ind_2_sex == "M", "MomSon", # Mother-son
                                                       "UR")) # Unrelated
  
  
  MOP <- pop.pairwise.df_filt_female %>%
    filter(KinPair == "MomDau" | KinPair == "MomSon") %>%
    dplyr::select(Ind_1, Ind_1_sex, Ind_2, Ind_2_sex, KinPair)
  
  table(MOP$KinPair)
  
  
  # Retain minimum number of columns for the sample dataframe
  samples <- samples %>%
    dplyr::select(ID, birth.year, sampling.age, sex, lon, lat, sampling.year)
  
  
  # Save sample and MOPs
  save(samples, MOP, file = "michigan_fake_samples.RData")
  
 
  ### GET INDICES OF RELATIVE ABUNDANCE 
  
  # Look at your pop
  pop <- loopy.list[[estimation_year]]

  # Load Michigan counties in UP (broad RAI)
  UP_counties <- st_read("C:/Users/aseveque/Senckenberg Dropbox/Anthony Seveque/Postdoc_MSU/3. michigan_bear_application/Data/GIS/Counties_(v17a)/Counties_(v17a).shp") %>% # Upper Peninsula borders
    st_transform(32616) %>%
    filter(PENINSULA == "upper")
  
  # Sum number of bears in each county
  intersections <- st_intersects(UP_counties, pop)
  UP_counties$count <- lengths(intersections)
  
  # Also create a perfect RAI, just for proof of concept
  intersections <- st_intersects(sf_grid_small, pop)
  sf_grid_small$count <- lengths(intersections)
  
  # Get relative abundance
  density_mean <- mean(UP_counties$count) # Get the mean of population density in your study area
  UP_counties$rel_count <- (UP_counties$count) / density_mean # Cross-multiplication (rule of three), so that basically new values as if mean = 1
  
  density_mean <- mean(sf_grid_small$count)
  sf_grid_small$rel_count <- (sf_grid_small$count) / density_mean
  
  # Rasterize the multipolygon
  template_raster_UTM <- rast(xmin= 238021.9, ymin= 4993636, xmax= 740236.1, ymax= 5259001, crs = "epsg:32616")
  UP_counties_raster <- rasterize(UP_counties, template_raster_UTM, field = "rel_count")
  
  template_raster_UTM <- rast(xmin= 238021.9, ymin= 4993636, xmax= 740236.1, ymax= 5259001, crs = "epsg:32616")
  sf_grid_small_rasterized <- rasterize(sf_grid_small, template_raster_UTM, field = "rel_count")
  
  # Save
  
  broad_raster_RDI <- UP_counties_raster
  fine_raster_RDI <- sf_grid_small_rasterized
  
  plot(broad_raster_RDI)
  plot(fine_raster_RDI)

  writeRaster(broad_raster_RAI, "michigan_broad_RDI.tif", overwrite=TRUE)
  writeRaster(fine_raster_RAI, "michigan_fine_RDI.tif", overwrite=TRUE)
  
  # POPULATION KDE
  
  # Convert sf points to a data frame with coordinates
  point_coords <- st_coordinates(pop)
  df <- data.frame(lon = point_coords[,1], lat = point_coords[,2])

  # Create the density plot
  # ggplot(df, aes(x = lon, y = lat)) +
  #   stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.7) +
  #   geom_point(color = "red", alpha = 0.5, size = 1) +
  #   scale_fill_viridis() +
  #   theme_minimal() +
  #   labs(x = "Longitude", y = "Latitude", title = "pop", fill = "Density")  
  
  # SAMPLE KDE
  # Convert sf points to a data frame with coordinates
  point_coords <- st_coordinates(samples.sf)
  df <- data.frame(lon = point_coords[,1], lat = point_coords[,2])
  
  # Create the density plot
  # ggplot(df, aes(x = lon, y = lat)) +
  #   stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.7) +
  #   geom_point(color = "red", alpha = 0.5, size = 1) +
  #   scale_fill_viridis() +
  #   theme_minimal() +
  #   labs(x = "Longitude", y = "Latitude", title = "samples", fill = "Density")
  