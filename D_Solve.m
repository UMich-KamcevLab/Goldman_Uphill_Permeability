function [D_fit, D_conf]=D_Solve
%Iterates over numerous simulations of a permeability experiment to regress
%for the diffusivities of each ion which best capture the experimental data
%   Specify the input file, which should be formatted in the template of
%   the provided Permeability_Data.xlsx for compatability with
%   Permeability_Import.m
%
%   Output will be diffusivities for Li, Mg, and Cl in cm^2/s, alongside a
%       matrix of the lower and upper bounds of a 95% confidence interval
%       around the regressed values, calculated based on the convergence of
%       the nonlinear solver.

[filename, path]=uigetfile('*.xlsx');
file=fullfile(path,filename);

%Import data from spreasheet and unpack experimental dataset
[experiment, cells, membrane] = Permeability_Import(file);
time_expt = experiment.Time; %min
C_Li_expt = experiment.ReceiverLithium; %kM
C_Mg_expt = experiment.ReceiverMagneisum; %kM

%Calculate the duration of the simulation needed to encompass the data
time_sim = round(time_expt(end)+30,-1);

%Establish the bounds and initial guess for ion diffusivities
D_min = zeros(3,1);
D_max = D_min + 1;
D_initial = D_min + 0.1;

%Set up options for the nonlinear solver
%Algorithm selection
alg = 'levenberg-marquardt';
%alg = 'trust-region-reflective';

%Select tolerance and stop criteria
options = optimoptions('lsqnonlin', 'Algorithm', alg, 'Display', 'iter',...
    'MaxIter', 500, 'MaxFunctionEvaluations', 10^3, ...
    'TolFun', 10^-4, 'Tolx', 10^-8);

%Solve for ion diffusivities
%Bounds are inclusive, so ±10^-10 makes the bounds exclusive
%output the fited values, D_fit, the residuals, R, and the jacobian, J
[D_fit, ~, R, ~, ~, ~, J]= lsqnonlin(@Sim_Error, D_initial, D_min+10^-10, D_max-10^-10, options);

%Calculate the 95% confidence interval around the fitted results
D_conf=nlparci(D_fit, R, 'Jacobian', J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot comparisons of the simulated results and the experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate the concentration profile of the final result
C_rec = Perm_Sim(D_fit, time_sim, membrane, cells, [0,0]);

figure;
subplot(2,1,1);
title(membrane.Name +" "+ cells.Configuration)
hold on
plot(C_rec{:,'time'},C_rec{:,"Li"},'b-')
plot(C_rec{:,'time'},C_rec{:,"Mg"},'r-')
plot(time_expt,C_Li_expt,'bx')
plot(time_expt,C_Mg_expt,'rx')
hold off
xlabel("Time (min)")
ylabel("Concentration (mM)")
legend("Li^{+} Simulation", "Mg^{2+} Simulation","Li^{+} Experiment", "Mg^{2+} Experiment",'location','northwest')

subplot(2,1,2);
hold on
plot(C_rec{:,'time'},C_rec{:,"Li"}./C_rec{:,"Mg"},'-')
plot(time_expt,C_Li_expt./C_Mg_expt,'x')
hold off
ylim([1,max(C_Li_expt./C_Mg_expt)*1.2])
xlabel("Time (min)")
ylabel("Separation Factor")
legend("Simulation", "Experiment", 'location','northeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prepare diffusivities for output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert dimensionless diffusivities to cm^2/s
[D_fit(1),D_fit(2),D_fit(3)] = DiffusivityDimensions(D_fit);
[D_conf(1,1),D_conf(2,1),D_conf(3,1)] = DiffusivityDimensions(D_conf(:,1));
[D_conf(1,2),D_conf(2,2),D_conf(3,2)] = DiffusivityDimensions(D_conf(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Error function which is set to be minimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Error = Sim_Error(D_ratio)
    %Simulates permeabilities given the vector D_ratio
    %Outputs a vector of normalized residuals for the concentrations

    %Run Simulation using D_ratio
    C_rec = Perm_Sim(D_ratio, time_sim, membrane, cells, [0,0]);

    %Interpolate results to the experimental times
    C_Li_interp = interp1(C_rec{:,'time'}, C_rec{:,'Li'}, time_expt);
    C_Mg_interp = interp1(C_rec{:,'time'}, C_rec{:,'Mg'}, time_expt);
    
    %Calculate the residual for each ion, weighting Mg by half an order of
    %magnitude to account for the lower concentrations
    NormResid_Li = (C_Li_expt - C_Li_interp);
    NormResid_Mg = (C_Mg_expt - C_Mg_interp)*3.33;
    Error = [NormResid_Li; NormResid_Mg];
end
end
