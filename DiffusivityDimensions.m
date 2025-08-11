function [D_Li_mem, D_Mg_mem, D_Cl_mem] = DiffusivityDimensions(D_ratio)
%Converts normalized diffusivities for Li, Mg, Cl to membrane values
%   Input a dimensionless vector of normalized diffusivites for Li, Mg, Cl
%       This ratio represents D_mem / D_sol for Li and Cl
%       Divalents are broadly slower than monovalents in membranes. To gain
%           additional resolution, the value for Mg is normalized against
%           the normalized Lithium diffusivity too.
%   The function will output membrane diffusivities in cm^2/s
%   The calculation uses solution diffusivities calculated from the
%       limiting equivalent conductances reported in the Robinson & Stokes 
%       textbook, 1965, at 25 ºC

%Define constants and experimetnal temperature
F_const = 96500; %Faraday's constant (C/mol)
R_const = 8.314; %Gas constant (J/(mol*K))
T = 273.15+25; %Absolute temperature (K)

%Define ion valences
z_Li = 1; % Valence of Li+
z_Mg = 2; % Valence of Mg2+
z_Cl = -1; %Valence of Cl-

%Convert limiting equivalent conductances to diffusivities using the
%Nernst-Einstein equation, D = Lambda_0 * R*T/z/F^2
D_Li_sol = 38.68/abs(z_Li)/F_const^2*R_const*T; %cm^2/s
D_Mg_sol = 53.05/abs(z_Mg)/F_const^2*R_const*T; %cm^2/s
D_Cl_sol = 76.35/abs(z_Cl)/F_const^2*R_const*T; %cm^2/s

%Membrane
%Set norm diffusivity of magnesium to max at that of lithium
D_ratio(2)=D_ratio(1)*D_ratio(2);

D_Li_mem = D_ratio(1)*D_Li_sol; %cm^2/s
D_Mg_mem = D_ratio(2)*D_Mg_sol; %cm^2/s
D_Cl_mem = D_ratio(3)*D_Cl_sol; %cm^2/s
end