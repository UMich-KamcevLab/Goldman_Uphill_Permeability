function [C_rec] = Perm_Sim(D_ratio, time, membrane, cells, genplot)
%Simulates LiCl/MgCl mixed salt permeability experiment for specified
%normalized diffusivities of each ion, outputting the downstream
%concentration evolution over time
%   Function input is a dimensionless membrane diffusivity for each of
%       Li+, Mg2+, and Cl-, in that order. Values are dimensionalized using
%       the function DiffusivityDimensions.m.
%   Additionally, provide the time duration in minutes for the duration of
%       the simulation. 
%   Finally, provide structs containing membrane and cell
%       configurations, including concentrations, thicknesses, areas, and
%       volumes. These structs are acquired from Permeability_Import.m
%
%   Function output is the concentrations of Li+, Mg2+, and Cl- (mM) as a
%       function of time (min) in the downstream receiver chamber. Output
%       is formatted as a table. 
%
%   To handle startup, because outputs are not very reliable for the first
%       few simulated timepoints, they are excluded from the output. 
%
%   Donnan equiulibrium is handled via a separate function specifically
%       tuned to handle Li+, Mg2+, Cl- equilibrium, assuming the Ideal
%       Donnan model holds at the downstream side of the membrane:
%       Donnan_LiMgCl.m
%
%   The input genplot toggles the production of graphical results. It
%       should be a 2 point matrix of 0 (off) or 1 (on).
%   Genplot Key
%       The first bool value is the plot of concentrations over time
%       The second bool value is the plot of fluxes and SF over time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import Settings and Establish Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ion Valences
z_Li = 1; % Valence of Li+
z_Mg = 2; % Valence of Mg2+
z_Cl = -1; %Valence of Cl-

%Simulation Settings
time_max=time*60; %convert min to s
dt=10; %resoltuion of simulation, s

%%%%%%%%%%
%Preallocate times, concentrations and fluxes
%Time vector
time_steps = (0:dt:time_max).'; %Vector of simulation times, s
num_steps=length(time_steps); %Number of simulated steps

%Concentration vectors
C_rec_Li = nan(num_steps, 1); % Receiver Li+ concentration over time
C_rec_Mg = nan(num_steps, 1); % Receiver Mg2+ concentration over time
C_rec_Cl = nan(num_steps,1); %Receiver Cl- concentration over time
C_mem_Li_rec = nan(num_steps, 1); % Membrane Li+ concentration over time (receiver side)
C_mem_Mg_rec = nan(num_steps, 1); % Membrane Mg2+ concentration over time (receiver side)
C_mem_Cl_rec = nan(num_steps, 1); % Membrane Cl- concentration over time (receiver side)

%Flux vectors
J_Li = nan(num_steps, 1); % Li+ flux over time
J_Mg = nan(num_steps, 1); % Mg2+ flux over time
J_Cl = nan(num_steps,1); %Cl in flux over time

%%%%%%%%%%
%Load in experimental parameters the functions
%Experiment Constants
V = cells.Volume; %cm^3
S = cells.SurfaceArea; %cm^2
C_rec_Li(1) = cells.ReceiverLithium; %kM
C_rec_Mg(1) = cells.ReceiverMagnesium; %kM

%Membrane Parameters
WVF = membrane.WaterVolumeFraction; %dimensionless
C_x = membrane.FixedChargeConcentration; %kM
d = membrane.Thickness; %cm
C_mem_Li_don = membrane.DonorLithium;
C_mem_Mg_don = membrane.DonorMagnesium;

%Membrane diffusivities
[D_Li, D_Mg, D_Cl]=DiffusivityDimensions(D_ratio);

%%%%%%%%%%
%Process the parameters if needed
%Can't work with starting downstream of DI water, because of Donnan solver
%Set these values to be very small if DI water is used downstream
if C_rec_Li(1)==0 && C_rec_Mg(1)==0
    C_rec_Li(1)=10^-10;
    C_rec_Mg(1)=10^-10;
end

%Set chloride concentrations by electroneutrality
C_rec_Cl(1) = (z_Li*C_rec_Li(1)+z_Mg*C_rec_Mg(1))/-z_Cl;
C_mem_Cl_don = (z_Li*C_mem_Li_don+z_Mg*C_mem_Mg_don-C_x)/-z_Cl;

%%%%%%%%%%
%Define functions to simulate fluxes and mass balancing

%Flux sub-functions from citation
A1 =@(C_mem_Li_rec) D_Li * C_mem_Li_rec + D_Cl * C_mem_Cl_don;
A2 =@(C_mem_Mg_rec) D_Mg * C_mem_Mg_rec;
B1 =@(C_mem_Cl_rec) D_Li * C_mem_Li_don + D_Cl * C_mem_Cl_rec;
B2 = D_Mg * C_mem_Mg_don;
P =@(C_mem_Li_rec, C_mem_Mg_rec, C_mem_Cl_rec) (A1(C_mem_Li_rec) - B1(C_mem_Cl_rec) + sqrt((B1(C_mem_Cl_rec) - A1(C_mem_Li_rec))^2 + 4*(B1(C_mem_Cl_rec) + 4*B2)*(A1(C_mem_Li_rec) + 4*A2(C_mem_Mg_rec))))/(2*(B1(C_mem_Cl_rec)+4*B2));

%Use P to calculate flux, J, mol/cm^2/s
J =@(P, z, D, C_mem_don, C_mem_rec) (z*D/d) * log(P) * -WVF*(C_mem_rec - C_mem_don * P^z)/(P^z-1);

%Use Flux to update receiver concentration
C_rec_new =@(C_rec_old, J_old) C_rec_old + J_old *S/V*dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Execute Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Perform initialization of permeability iteration
t=1;

%Establish Equilibrium with initial solution
[C_mem_Li_rec(t), C_mem_Mg_rec(t), C_mem_Cl_rec(t)] = Donnan_LiMgCl(C_rec_Li(t), C_rec_Mg(t), C_rec_Cl(t), C_x);

%Determine initial fluxes
P_current=P(C_mem_Li_rec(t), C_mem_Mg_rec(t), C_mem_Cl_rec(t));
J_Li(t)=J(P_current, z_Li, D_Li, C_mem_Li_don, C_mem_Li_rec(t));
J_Mg(t)=J(P_current, z_Mg, D_Mg, C_mem_Mg_don, C_mem_Mg_rec(t));
J_Cl(t)=J(P_current, z_Cl, D_Cl, C_mem_Cl_don, C_mem_Cl_rec(t));

%Iterate for all times
for t=2:num_steps

    %Update receiver cell concentration
    %Use previous time flux
    C_rec_Li(t)=max(C_rec_new(C_rec_Li(t-1), J_Li(t-1)), 0);
    C_rec_Mg(t)=max(C_rec_new(C_rec_Mg(t-1), J_Mg(t-1)), 0);
    C_rec_Cl(t)=max(C_rec_new(C_rec_Cl(t-1), J_Cl(t-1)), 0);
    %Max functions ensure that, if approaching 0, 
    %concentrations don't drop below 0.


    %Update receiver-side membrane concentration
    %Establish Donnan equilibrium with new cell concentration
    [C_mem_Li_rec(t), C_mem_Mg_rec(t), C_mem_Cl_rec(t)] = Donnan_LiMgCl(C_rec_Li(t), C_rec_Mg(t), C_rec_Cl(t), C_x);

    %Calculate new fluxes
    %Use updated membrane concentrations
    P_current=P(C_mem_Li_rec(t), C_mem_Mg_rec(t), C_mem_Cl_rec(t));
    J_Li(t)=J(P_current, z_Li, D_Li, C_mem_Li_don, C_mem_Li_rec(t));
    J_Mg(t)=J(P_current, z_Mg, D_Mg, C_mem_Mg_don, C_mem_Mg_rec(t));
    J_Cl(t)=J(P_current, z_Cl, D_Cl, C_mem_Cl_don, C_mem_Cl_rec(t));
end

%Separation factor is ratio of lithium to magnesium
sep_factor=C_rec_Li(2:end)./C_rec_Mg(2:end);
sep_factor=[1; sep_factor];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Adjust Units and Prepare table for export
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_steps = time_steps ./60; %convert s to min
C_rec_Li = C_rec_Li .*10^6; %convert kM to mM
C_rec_Mg = C_rec_Mg .*10^6; %convert kM to mM
C_rec_Cl = C_rec_Cl .*10^6; %convert kM to mM

%Export ignores the first 2 time steps
C_rec = table(time_steps(3:end), C_rec_Li(3:end), C_rec_Mg(3:end), C_rec_Cl(3:end));
C_rec.Properties.VariableNames=["time", "Li", "Mg", "Cl"];
C_rec.Properties.VariableUnits=["min","mM","mM","mM"];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate plots of the receiver cell, fluxes, and separation factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Concentrations of Receiver Cell
if genplot(1)
    figure;
    subplot(3,1,1);
    plot(C_rec{:,"time"}, C_rec{:,"Li"}, '-o');
    xlabel('Time (min)');
    ylabel('Li+ Conc. (mM)');
    title('Li+ Concentration in Receiver Chamber');
    
    subplot(3,1,2);
    plot(C_rec{:,"time"}, C_rec{:,"Mg"}, '-o');
    xlabel('Time (min)');
    ylabel('Mg2+ Conc. (mM)');
    title('Mg2+ Concentration in Receiver Chamber');
    
    subplot(3,1,3);
    plot(C_rec{:,"time"}, C_rec{:,"Cl"}, '-o');
    xlabel('Time (min)');
    ylabel('Cl- Conc. (mM)');
    title('Cl- Concentration in Receiver Chamber');
end

%Plot Fluxes and separation factor for Li+ / Mg2+
if genplot(2)
    figure
    subplot(2,1,1);
    hold on
    plot(C_rec{:,"time"}, J_Li(3:end), '-o');
    plot(C_rec{:,"time"}, J_Mg(3:end), '-o');
    plot(C_rec{:,"time"}, J_Cl(3:end), '-o');
    hold off
    xlabel('Time (min)');
    ylabel('Flux (mol/cm^{2}/s)');
    legend('Li^{+}', 'Mg^{2+}', 'Cl^{–}')
    title('Ion Fluxes')
    
    subplot(2,1,2);
    plot(C_rec{:,"time"}, sep_factor(3:end), '-o');
    xlabel('Time (min)');
    ylabel('Separation Factor');
    title('Separation Factor');
end
end