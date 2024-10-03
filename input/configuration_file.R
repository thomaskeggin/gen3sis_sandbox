# General settings -------------------------------------------------------------
# set the random seed for the simulation
random_seed = 42

# set the starting time step or leave NA to use the earliest/highest timestep
start_time = NA

# set the end time step or leave as NA to use the lates/lowest timestep (0)
end_time = NA

# maximum total number of species in the simulation before it is aborted
max_number_of_species = 1000

# maximum number of species within one cell before the simulation is aborted
max_number_of_coexisting_species = 100

# a list of traits to include with each species
trait_names = c("thermal_optimum",
                "thermal_standard_deviation",
                "depth_limit")

# ranges to scale the input environments with:
# not listed variable:         no scaling takes place
# listed, set to NA:           the environmental variable will be scaled from [min, max] to [0, 1]
# listed with a given range r: the environmental variable will be scaled from [r1, r2] to [0, 1]
environmental_ranges = list()

# a place to inspect the internal state of the simulation and collect additional information if desired
end_of_timestep_observer = function(data, vars, config){
  save_richness()
  save_species()
  save_phylogeny()
}

# Initialisation ---------------------------------------------------------------
# the initial abundance of a newly colonized cell, both during setup and later when colonizing a cell during the dispersal
initial_abundance = 0.1

# place species within simulation:
create_ancestor_species <- function(landscape, config) {
  
  #browser()
  
  seascape <- 
    landscape$environment |> 
    as_tibble(rownames = "cell")
  
  coords <-
    landscape$coordinates |> 
    as_tibble(rownames = "cell")
  
  seascape <-
    left_join(seascape,coords,
              by = "cell")
  
  # define the starting cells for the two species --------
  # species starting in the Atlantic, limited by depth
  Pa_start_cells <-
    seascape |>
    filter(x > -88,x < -84, y > 20, depth > -1000) |>
    pull(cell)
  
  # species starting in the Pacific, not limited by depth
  Pp_start_cells <-
    seascape |>
    filter(x > -88,x < -84, y < 10) |>
    pull(cell)
  
  # Remember, the species object is just a list!
  species_object <- list()
  
  # create the Atlantic species-----------------------------
  species_object[[1]] <-
    gen3sis::create_species(initial_cells = Pa_start_cells,
                            config = config)
  
  # generate mean thermal niche and standard deviation
  species_object[[1]]$traits[,"thermal_optimum"] <-
    18
  
  species_object[[1]]$traits[,"thermal_standard_deviation"] <-
    0.3
  
  # depth trait
  species_object[[1]]$traits[ , "depth_limit"]   <-
    -200
  
  # tag on a species name
  species_object[[1]]$lineage <-
    "Pisces_atlanticus"
  
  # create Pacific species ----------------------------------
  species_object[[2]] <-
    gen3sis::create_species(initial_cells = Pp_start_cells,
                            config = config)
  
  # generate mean thermal niche and standard deviation
  species_object[[2]]$traits[,"thermal_optimum"] <-
    21
  
  species_object[[2]]$traits[,"thermal_standard_deviation"] <-
    0.1
  
  # depth trait
  species_object[[2]]$traits[ , "depth_limit"]   <-
    -10000
  
  # tag on a species name
  species_object[[2]]$lineage <-
    "Pisces_pacificus"
  
  # output species object
  return(species_object)
  
}

# Speciation -------------------------------------------------------------------
# threshold for genetic distance after which a speciation event takes place.
# speciation after every timestep : 0.9.
# we are removing the speciation dynamic by setting the threshold to infinity.
divergence_threshold = 3

# factor by which the genetic distance is increased between geographically isolated population of a species
# can also be a matrix between the different population clusters
get_divergence_factor <- function(species, cluster_indices, landscape, config) {
  
  return(1)
}

# Dispersal --------------------------------------------------------------------
# returns n dispersal values
get_dispersal_values <- function(num_draws, species, landscape, config) {
  
  return(
    rweibull(num_draws,
             shape = 2,
             scale = 500)
  )
}

# Evolution --------------------------------------------------------------------
# mutate the traits of a species and return the new traits matrix
apply_evolution <- function(species, cluster_indices, landscape, config){
  
  #browser()
  
  traits <-
    species[["traits"]]
  
  traits[,"thermal_optimum"] <-
    traits[,"thermal_optimum"] + rnorm(length(traits[,"thermal_optimum"]),
                                       mean=0,
                                       sd=0.04)
  
  return(traits)
}

# Ecology ----------------------------------------------------------------------
# called for every cell with all occuring species, this function calculates who survives in the current cells
# returns a vector of abundances
# set the abundance to 0 for every species supposed to die
apply_ecology <- function(abundance, traits, local_environment, config) {
  
  #browser()
  
  new_abundance <-
    
    # trait information
    dplyr::as_tibble(traits) |> 
    
    # environmental information
    cbind(dplyr::as_tibble(local_environment)) |> 
    
    # ecology calculations
    dplyr::mutate(
      
      # start abundance
      start_abundance = abundance, 
      
      # the distribution density if the species' niche perfectly fits the temperature
      optimal_density = dnorm(thermal_optimum,
                              mean = thermal_optimum,
                              sd = thermal_standard_deviation),
      
      # the distribution density given the distance between species niche and temperature
      species_density = dnorm(thermal_optimum,
                              mean = temp,
                              sd = thermal_standard_deviation),
      
      # determine the abundance
      end_abundance = species_density/optimal_density,
      
      # drive to (local) extinction if abundance is below 10%
      end_abundance = ifelse(end_abundance < 0.1,
                             0,
                             end_abundance),
      
      # drive to (local) extinction if the temperature is completely unsuitable
      end_abundance = ifelse(species_density == 0,
                             0,
                             end_abundance),
      
      # drive to (local) extinction if the depth is unsuitable
      end_abundance = ifelse(depth < depth_limit,
                             0,
                             end_abundance)
    ) |> 
    
    # extract end abundance only
    dplyr::pull(end_abundance)
  
  # assign cell names
  names(new_abundance) <- names(abundance)
  
  # fin
  return(new_abundance)
}