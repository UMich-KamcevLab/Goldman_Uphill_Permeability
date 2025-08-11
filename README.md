# Goldman_Uphill_Permeability
This package of MATLAB functions simulates LiCl and MgCl2 permeation across ion-exchange membranes. Modeling is achieved in MATLAB by combining the Ideal Donnan model and the Goldman Equation into a time-resolved simulation. 

## System Requirements / Dependencies
In addition to a basic MATLAB license, the non-linear fitting requires the following toolboxes:  
- Optimization Toolbox  
- Statistics and Machine Learning Toolbox  

The code was developed in MATLAB R2023b.  

## Quick Use Case
Download all of the files to a single directory accessible via MATLAB.  
Run D_Solve.m and select the provided file to analyze, with no other inputs. The function will generate concentration profiles, fluxes, separation factors, and a fit of the experimental data.  
The output of the function is the diffusivity fit for each of Li+, Mg2+, and Cl- in the membrane in cm2/s.  
Depending on the duration of the epxeirment being simulated, D_Solve.m may take multiple minutes to execute.  

## Related Publications and Citation
This simulation platform was originally implemented by Higa *et al., Journal of Membrane Science*, 1988, 37, 3, 251-266, DOI: 10.1016/S0376-7388(00)82432-1.  
This code was developed for publication alongside Santiago-Pagán and Patel *et al., Nature Chemical Engineering*, 2025.  
To cite this package, please reference Santiago-Pagán and Patel *et al.*  

## Full Contents and Usage
Included are scripts to:  
- Import data from the excel template: Permeability_Import.m  
- Convert dimensionless diffusivities to values in cm^2/s: DiffusivityDimensions.m  
- Calculate interfacial equilibrium using the Ideal Donnan model: Donnan_LiMgCl.m  
- Simulate a permeability experiment for given experimental conditions: Perm_Sim.m  
- Non-linearly solve for the apparent diffusivities of each ion in the experiment: D_Solve.m  

Other files include:  
- A formatted data template to allow for easy import of experimental permeability data and accompanying cell / membrane information: Permeability_Data.xlsx  
- A DEMO folder containing example data that can be fit using this simulation: Permeability_Data_CR61Example.xlsx  