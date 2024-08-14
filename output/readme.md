### Output

This folder contains all the outputs from the model. In this case, we have just one simulation, but we could have many, so we have a directory for the outputs of each simulation. These are flexible, but the core outputs include:

- **Species**: directory containing the species objects (.rds) for each time step. Each one is a list of species, containing:
  - **id** - species ID
  - **abundance** - vector of abundances in occupied cells
  - **traits** - data frame of traits for each occupied cell
  - **divergence** - list containing two compression objects for the cell-to-cell divergence values between occupied cells
- **Configuration.R**: a copy of the configuration file for that simulation.
- **Starting_ranges.pdf**: plot of the species ranges at the start of the simulation.
- **Starting_richness.pdf**: plot of the species richness at the start of the simulation.