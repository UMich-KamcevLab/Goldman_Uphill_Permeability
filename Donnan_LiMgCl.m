function [C_mem_Li, C_mem_Mg, C_mem_Cl] = Donnan_LiMgCl(C_sol_Li, C_sol_Mg, C_sol_Cl, C_x)
%Solves Ideal Donnan equilibrium for the given conditions
%Specifically tuned for CEMs contacting mixtures of LiCl and MgCl2 salts
%   Input external solution concentrations as kM (kmol/L [=] mol/mL)
%   Input membrane fixed charge concentration (FCC, Cx) as kM
%   Input concentration should be mol fixed charge / mL sorbed H2O
%   Output concentrations will also be normalized to mL sorbed H2O
%   Solving method will isolate the Donnan potential from the
%       characteristic polynomial, of degree max(z_cation)-min(z_anion) = 3

%Solve for the exponential of the Donnan Potential (Psi) nondimensionalized by 
% the Faraday constant (F), ideal gas constant (R), and absolute temperature (T):
%exp(-F/R/T*Psi), which must be greater than 0.
Donnan_exp = roots([2*C_sol_Mg, C_sol_Li, -C_x, -C_sol_Cl]); %solve the cubic
Donnan_exp = Donnan_exp(Donnan_exp>0); %select the valid result

if length(Donnan_exp)>1
    error("Unable to discern between multiple Donnan potential possibilities.\n")
end

%Solve interfacial equilibrium for each species
C_mem_Li = C_sol_Li*Donnan_exp;
C_mem_Mg = C_sol_Mg*Donnan_exp^2;
C_mem_Cl = C_sol_Cl/Donnan_exp;
end