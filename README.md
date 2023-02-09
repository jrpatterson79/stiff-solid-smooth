## Overview
The code in this repository is used to generate the analysis presented in Patterson and Cardiff (2023), _Stiff, Solid, and Smooth?: Complex Fracture Hydraulics Revealed through Oscillatory Flow Interference Testing_, a manuscript submitted to Water Resources Research for peer review. This paper describes a suite of numerical oscillatory hydraulic tomography experiments to explore multiple mechanisms that lead to period-dependent fracture flow parameters when characterizing bedrock fractures using oscillatory pressure signals. 

## Models
This folder contains MATLAB scripts for each model used to simulate oscillatory flow tests and generate synthetic data used for parameter estimation. 
* fracture_heterogeneity.m: 2-D numerical model that explores the impact of fracture heterogeneity on parameter estimates
* validation_2d.m: Uses analytical model to validate the above numerical model
* fracture_host_rock_fluid_exchange.m: 3-D numerical model that explores the impact of fluid exchange between the fracture and host rock on parameter estimates.
* validation_3d.m: Uses analytical model to validate the above analytical model
* fracture_hydromechanics.mph: COMSOL hydro mechanical model to generate pressure and fracture displacement signals and output data as .csv file
* hydromechanics_pest.m: Uses COMSOL outputs to conduct parameter estimation 
All MATLAB scripts were developed using MATLAB 2019b and have recently been executed using MATLAB 2021b. Hydromechanical modeling was conducted using COMSOL 5.6

## func_lib
This folder contains all of the necessary function files to execute the modeling scripts described above.

## figure_scripts
This folder contains MATLAB scripts that reproduce figures presented throughout the manuscript.

## files
This folder contains the input files needed to execute the scripts in *models* and *figure_scripts*

## License
The code and data are provided as open source under the GNU General Public License v3.0. It is provided without warranty, but should perform as described in the manuscript when executed without modification. If you would like to use any of the code in this repository for research, software, or publications, we ask that you provide a citation to the code and journal article (See references below).