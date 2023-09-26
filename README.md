## Overview
The code in this repository is used to generate the analysis presented in Patterson and Cardiff (2023), _Stiff, Smooth, and Solid?: Complex Fracture Hydraulics' Imprint on Oscillatory Flow Interference Testing_, a manuscript submitted to Water Resources Research for peer review. This paper describes a suite of numerical oscillatory hydraulic tomography (OHT) experiments to explore potential mechanisms contributing to an period-dependence in effective fracture hydraulic parameter estimates returned when characterizing bedrock fractures using oscillatory pressure signals. 

## models
This folder contains MATLAB scripts for each model used to simulate oscillatory flow tests and generate synthetic data used for parameter estimation. 
* fracture_heterogeneity.m: 2-D numerical groundwater flow model that explores the impact of aperture heterogeneity on effective fracture hydraulic parameter estimates
* validation_2d.m: This model verifies the accuracy of the above numerical model by comparing head phasors simulated by a homogeneous 2-D numerical model to head phasors simulated for a fully confined aquifer as described by Rasmussen et al. (2003)
* fracture_host_rock_fluid_exchange.m: 3-D numerical groundwater flow model that explores the impact of fluid exchange between the fracture and host rock on parameter estimates
* validation_3d.m: This model verifies the accuracy of the above numerical model by comparing head phasors simulated by a homogeneous 3-D numerical model to head phasors simulated for a leaky confined aquifer as described by Rasmussen et al. (2003)
* fracture_hydromechanics.mph: Axi-symmeteric COMSOL hydromechanical model that simulates pressure and fracture displacement signals and outputs the simulated data as .csv file, which serves as the input for hydromechanics_pest.m
* hydromechanics_pest.m: Uses COMSOL outputs to conduct effective fracture hydraulic parameter estimation
 
_All MATLAB scripts were developed using MATLAB 2019b and have recently been executed using MATLAB 2021b. Hydromechanical modeling was conducted using COMSOL 5.6_

## func_lib
This folder contains all of the necessary function files to execute the modeling scripts described above.

## figure_scripts
This folder contains MATLAB scripts that reproduce figures presented throughout the manuscript.

## files
This folder contains the input files needed to execute the MATLAB scripts in *models* and *figure_scripts*

## License
The code and data are provided as open source under the GNU General Public License v3.0. It is provided without warranty, but should perform as described in the manuscript when executed without modification. If you would like to use any of the code in this repository for research, software, or publications, we ask that you provide a citation to the code and journal article (See references below).

## References
Patterson, J.R., Cardiff, M. (2023). Stiff, Smooth, and Solid?: Complex Fracture Hydraulics' Imprint on Oscillatory Flow Interference Testing. *Under review for Water Resources Research*

Rasmussen, T. C., Haborak, K. G., & Young, M. H. (2003). Estimating aquifer hydraulic properties using sinusoidal pumping at the Savannah River site, South Carolina, USA. Hydrogeology Journal, 11(4), 466–482. https://doi.org/10.1007/s10040-003-0255-7
