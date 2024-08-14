### Input

We are using a very coarse resolution representation of shallow water marine reefs. The processing steps to get these data were complex and the external data quite large, so we are skipping over these steps. We start here with the pre-formatted landscape object and distance matrices.

This directory contains the model inputs:

- **landscape.rds**: this is a list of xy data frames containing the environmental variables for each time step.
- **distance_matrices_full**: the geographic distances between each habitable cell - computed from the landscape (seascape) object.