function [experiment, cells, membrane] = Permeability_Import(file)
%Imports experimental data of LiCl/MgCl2 permeability experiments
%   Follows the template of the included file Permeability_Data.xlsx
%       Extracts the experimental concentration profile of the downstream
%       chamber, the conditions of the experiment, and the properties of
%       the membrane that are needed to simulate the permeability
%       experiment. 
%   Converts the inputs from units that are easier to measure and input
%       into units easier to simulate, i.e., converting mM to kM[=]mol/cm^3

A = readtable(file,'NumHeaderLines',2);

%Labels
membrane.Name = A{1,13}{1}; %string
cells.Configuration = A{1,5}{1}; %string

%Sizes
cells.Volume = A{1,6}; %cm^3
cells.SurfaceArea = A{1,7}; %cm^2
membrane.Thickness = A{1,14}; %um
%
membrane.Thickness = membrane.Thickness/10000; %convert um to cm

%Membrane Properties
membrane.WaterVolumeFraction = A{1,15}; %dimensionless
membrane.FixedChargeConcentration = A{1,16}; %mol[fixed charge] / L[sorbed water]
membrane.DonorLithium = A{1,17}; %mol[counter-ion] / L[sorbed water]
membrane.DonorMagnesium = A{1,18}; %mol[counter-ion] / L[sorbed water]
%
membrane.FixedChargeConcentration = membrane.FixedChargeConcentration/1000; %convert M to kM
membrane.DonorLithium = membrane.DonorLithium/1000; %convert M to kM
membrane.DonorMagnesium = membrane.DonorMagnesium/1000; %convert M to kM

%Cell Concentrations
cells.DonorLithium = A{1,8}; %mmol / L[solution]
cells.DononrMagnesium = A{1,9}; %mmol / L[solution]
cells.ReceiverLithium = A{1,10}; %mmol / L[solution]
cells.ReceiverMagnesium = A{1,11}; %mmol / L[solution]
%
cells.DonorLithium = cells.DonorLithium/10^6; %convert mM to kM
cells.DononrMagnesium = cells.DononrMagnesium/10^6; %convert mM to kM
cells.ReceiverLithium = cells.ReceiverLithium/10^6; %convert mM to kM
cells.ReceiverMagnesium = cells.ReceiverMagnesium/10^6; %convert mM to kM

%Downstream Concentration Profile
experiment.Time = A{:,1}; %min
experiment.ReceiverLithium = A{:,2}; %mmol / L[solution]
experiment.ReceiverMagneisum = A{:,3}; %mmol / L[solution]
end