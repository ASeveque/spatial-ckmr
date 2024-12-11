####---------- Population density map ----------####

pop.density <- function(pop.sf){
  
  # Count how many individuals in each small grid cells
  # This is possible only because it is a simulation and we (the R user) are omniscient
  
  pop.sf <- pop.sf[pop.sf$age >= 1,] # look at density after dispersal
  
  # Create a box (100x100km)
  box <- matrix(c(0, 0, # Create the bounding points
                  0, 1,
                  1, 1,
                  1, 0,
                  0, 0),
                ncol=2, byrow=TRUE)
  
  sf_box <- st_polygon(list(box)) # Join the points to create a polygon
  sf_box <- st_sfc(sf_box, crs = 4326) # Transform polygon into a spatial feature
  
  # Create grid where animals are located
  sf_grid_small  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.05, 0.05)) %>% # Make it 50x50
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326)
  
  # Get centroid coordinates of each cell in the grid
  sf_grid_small_centroid <- st_centroid(sf_grid_small)
  
  # Count how many individuals are in each grid cell 
  # NB: this does NOT count the empty cells, need to add them manually
  pop.density.temp <- pop.sf %>%
    arrange(ID.small) %>%
    as.data.frame() %>%
    dplyr::select(ID.small) %>%
    plyr::count() 
  
  # Add empty cells with 0
  pop.density.temp2 <- data.frame(ID.small = c(1:400)) %>% # Empty cells are given an NA value
    merge(pop.density.temp, by="ID.small", all.x=TRUE)
  pop.density.temp2[is.na(pop.density.temp2)] <- 0 # Transform NAs into 0s
  
  # Create a raster with the relative densities (here, the absolute number of animals in each cells)
  # First, extract coordinates for each cell's centroid into a dataframe
  pop.xy <- st_coordinates (sf_grid_small_centroid)
  
  # Create 3D array with X (lon), Y (lat) and Z (counts)
  pop.density <- cbind (x = as.numeric(pop.xy[,1]),
                        y = as.numeric(pop.xy[,2]),
                        z = pop.density.temp2$freq)
  
  # Create a template raster just to have a crs and boundaries to align the data with. NOTE: this is based on the study area (sf_box).
  
  # First, as WGS 84
  template_raster <- rast(ncols = 20, nrows = 20, xmin= 0, ymin= 0, xmax= 1, ymax= 1, crs = "epsg:4326")
  pop.density.raster <- rasterize(pop.density[,1:2], template_raster, pop.density[,3])
  
  # Then, transform in UTM
  template_raster_UTM <- rast(ncols = 20, nrows = 20, xmin= 166021.4, ymin= 0, xmax= 277438.3, ymax= 110682.8, crs = "epsg:32631")
  pop.density.raster <- project (pop.density.raster, template_raster_UTM) # Change projection from 4326 (lon/Lat) to UTM)
  
  # Transform absolute value of abundance into relative abundance/density estimates
  density_mean <- as.numeric(global(pop.density.raster, "mean")) # Get the mean of population density in your study area
  density.pop.raster <- (1 * pop.density.raster) / density_mean # Cross-multiplication (rule of three), so that basically new values as if mean = 1
  
}


####---------- Samples density map ----------####

relative.density.samples <- function(samples){
  
  # Count how many individuals in each small grid cells
  # This is possible only because it is a simulation and we (the R user) are omniscient
  
  # Create a box (100x100km)
  box <- matrix(c(0, 0, # Create the bounding points
                  0, 1,
                  1, 1,
                  1, 0,
                  0, 0),
                ncol=2, byrow=TRUE)
  
  sf_box <- st_polygon(list(box)) # Join the points to create a polygon
  sf_box <- st_sfc(sf_box, crs = 4326) # Transform polygon into a spatial feature
  
  # Create grid where animals are located
  sf_grid_small  <- st_make_grid(sf_box, # Grid will follow the box's bounding box and crs (if any)
                                 square = TRUE,
                                 cellsize = c(0.05, 0.05)) %>% # Make it 50x50
    cbind(data.frame(ID = sprintf(paste("%0",nchar(length(.)),"d",sep=""), 1:length(.)))) %>% # Create a unique ID for each cell
    st_sf(crs = 4326)
  
  # Get centroid coordinates of each cell in the grid
  sf_grid_small_centroid <- st_centroid(sf_grid_small)
  
  samples$ID.small <- as.numeric(samples$ID.small)
  
  # Count how many sampled individuals are in each grid cell 
  # NB: this does NOT count the empty cells, need to add them manually
  sample.density.temp <- samples %>%
    arrange(ID.small) %>%
    as.data.frame() %>%
    dplyr::select(ID.small) %>%
    plyr::count()
  
  # Add empty cells with 0
  sample.density.temp2 <- data.frame(ID.small = c(1:400)) %>% # Empty cells are given an NA value
    merge(sample.density.temp, by="ID.small", all.x=TRUE)
  sample.density.temp2[is.na(sample.density.temp2)] <- 0 # Transform NAs into 0s
  
  # Create a raster with the relative densities (here, the absolute number of sampled animals in each cells)
  # First, extract coordinates for each cell's centroid into a dataframe
  sample.xy <- st_coordinates (sf_grid_small_centroid)
  
  # Create 3D array with X (lon), Y (lat) and Z (counts)
  sample.density <- cbind (x = as.numeric(sample.xy[,1]),
                           y = as.numeric(sample.xy[,2]),
                           z = sample.density.temp2$freq)
  
  # Create a template raster just to have a crs and boundaries to align the data with. NOTE: this is based on the study area (sf_box).
  # First, as WGS 84
  template_raster <- rast(ncols = 20, nrows = 20, xmin= 0, ymin= 0, xmax= 1, ymax= 1, crs = "epsg:4326")
  sample.density.raster <- rasterize(sample.density[,1:2], template_raster, sample.density[,3])
  
  # Then, transform in UTM
  template_raster_UTM <- rast(ncols = 20, nrows = 20, xmin= 166021.4, ymin= 0, xmax= 277438.3, ymax= 110682.8, crs = "epsg:32631")
  sample.density.raster <- project (sample.density.raster, template_raster_UTM) # Change projection from 4326 (lon/Lat) to UTM)
  
  # Transform absolute value of abundance into relative abundance/density estimates
  density_mean <- as.numeric(global(sample.density.raster, "mean")) # Get the mean of population density in your study area
  density.samples.raster <- (1 * sample.density.raster) / density_mean # Cross product (basically new values as if mean = 1)
  
}