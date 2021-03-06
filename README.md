## What is the BFM

The **Biogeochemical Flux Model (BFM)** is a numerical model for the simulation 
of the dynamics of major biogeochemical properties in marine ecosystems.

The BFM is open source software freely available under the GNU Public License. 

The model description is available on the BFM website http://www.bfm-community.eu and in the doc/ directory. Instructions for compilation and running are available on the web site http://www.bfm-community.eu/documentation

#### Members of the BFM agreement

The development of the BFM is supported by a system team through a formal agreement. The agreement establishes a program for model development that is discussed and prepared every year by the participants. Model developments are inserted by the system team into the different releases of the model. 
The agreement is open to entrance of more contributors and more information can be requested to info(at)bfm-community.eu
Funding members of the consortium are: CMCC, OGS and UNIBO.


## Biogeochemical Flux Model with the 1D Princeton Ocean Model

This repository contains the code of the **BFM_POM** coupled model to perform numerical experiments corresponding to the Sargasso Sea (BATS) by producing 10 years of daily data.

Two distinct BFM configurations with increasing complexity are made available.

**BFM56** is the full Biogeochemical Flux Model and includes:
- 4 phytoplankton groups
- 2 mesozooplankton groups
- 2 microzooplankton groups
- 1 bacteria group
- 5 nutrients
- particulate organic matter
- dissolved organic matter
- oxygen

**BFM17** is the reduced Biogeochemical Flux Model and includes:
- 1 phytoplankton group
- 1 microzooplankton group
- 3 nutrients
- particulate organic matter
- dissolved organic matter
- oxygen

### Environment setup
Both models have been set up to run using:
- gfortran
- openmpi
- netcdf

Please make sure that these have been correctly set up on your computer before proceeding.

### Configuration files

The parameter values for BFM56 are set in the namelists located in build/configurations/BFM56_POM1D. The parameter values for BFM17 are set in the namelists located in build/configurations/BFM17_POM1D. Before running either of the models, you may need to change the path assigned to the NETCDF variable in the run file. The path should correspond to the location of netCDF on your machine.

### Model execution

To run **BFM56+POM1D**, run the bash script in the main directory using:

```$> ./run_BFM56_POM.sh```

and output file with the model data will be stored in `./bfm_run/bfm56_pom1d/bfm56_pom1d.nc`

To run **BFM17+POM1D**, run the bash script in the main directory using:

```$> ./run_BFM17_POM.sh```

and output file with the model data will be stored in `./bfm_run/bfm17_pom1d/bfm17_pom1d.nc`
