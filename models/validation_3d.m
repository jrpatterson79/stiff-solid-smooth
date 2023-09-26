% Oscillatory Flow Testing Fracture Characterization
% Numerical Modeling Analysis

% This code validates the 3D numercial model used to explore the impacts of
% fracture-host rock fluid exchange on oscillatory flow signals. The numerical model is
% validated using the leaky confined analytical model developed by
% Rasmussen et al. (2003).

% Code developed by Jeremy Patterson
% Created Dec 2018; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Specify Directory
addpath(genpath('/.../.../')) % Provide location to func_lib subdirectory which contains needed function files

%% Forward Model Setup
% Specify Domain Geometry
xmin = -400; xmax = -xmin; dx = 2;
ymin = xmin; ymax = xmax; dy = dx;
dz = [5e-2; 0.075; 0.1125; 0.17; 0.25]; zmin = 0; zmax = sum(dz);

domain = struct('x',[],'y',[],'z',[]);
domain.x = [xmin : dx : xmax];
domain.y = [ymin : dy : ymax];
domain.z = [zmin: dz : zmax];

numx = numel(domain.x) - 1;
numy = numel(domain.y) - 1;
numz = numel(domain.z) - 1;
num_cells = numx * numy * numz;

% Grid Domain
[coords, cgrid] = plaid_cellcenter_coord(domain);
node_list = unique(coords(:,3));

f_idx = find(coords(:,3) == node_list(1));
c_idx = [find(coords(:,3) == node_list(2));...
         find(coords(:,3) == node_list(3))];

fx_coords = [coords(f_idx,1) coords(f_idx,2)];

% Specify boundary types and boundary values (x / y constant head
% boundaries, no flux in z)
bdry_types = [1; 1; 1; 1; 0; 0];
bdry_vals = zeros(6,1);
bdry_L = zeros(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals,'leaks',bdry_L);

% Specify Well Locations
elev = dz / 2; 
well_locs = [0 0 elev;...
             0 -4 elev;...
             -6 0 elev;...
             0 8 elev;...
            10 0 elev];
num_wells = numel(well_locs(:,1));

%% Generate test list
P = logspace(1, 3, 5);  % Oscillation period (s)
Q_max = 7e-5;           % Peak flow rate (m^3/s)

%Test list for OHT input [omega pump_well_index peak_flow_rate obs_well_index]
test_list = [];
r = [];
for i = 1:numel(P)
    for j = 2 : num_wells
        r  = [r; 
              sqrt((well_locs(1,1)-well_locs(j,1)).^2 + (well_locs(1,2)-well_locs(j,2)).^2)];
        test_list = [test_list; ...
                     (2*pi)/P(i) 1 Q_max j];
    end
end
num_tests = numel(test_list(:,1));

%% Specify Model Parameters
% Fluid and Rock Properties
rho = 998.23;           % Fluid density (kg/m^3)
g = 9.81;               % Gravitational acceleration (m/s^2)
mu = 1.0016e-3;         % Fluid dynamic viscosity (Pa s)
RockComp = 1.41e-8;     % Rock compressibility (Pa^-1)
WaterComp = 1 / 2.1e9;  % Fluid compressibility (Pa^-1)
eta_r = 0.15;           % Host rock porosity (-)

% Fracture Aperture Properties
ln_aper = log(3e-4); 
T_frac = ((rho*g)/(12 * mu)) * exp(ln_aper)^3;
K_frac = log(T_frac / (dz(1)*2));

Ss_rock = rho * g * (RockComp + (eta_r * WaterComp)); 
Ss_f = rho * g * WaterComp;
V_frac = Ss_f * dx * dy * 0.4;
V_rock = Ss_rock * dx * dy * 0.6;
Ss_eff = (V_frac + V_rock) / (dx*dy);
S_frac = Ss_eff * exp(ln_aper);
Ss_frac = log(S_frac / (dz(1)*2)); 

lnK = [K_frac log(1e-10) log(3e-1)]; %[Fracture Confining Bedrock]
lnSs =[Ss_frac -50 -1];              %[Fracture Confining Bedrock]

lnK_true = lnK(3) * ones(num_cells,1);   lnK_true(f_idx) = lnK(1);   lnK_true(c_idx) = lnK(2);
lnSs_true = lnSs(3) * ones(num_cells,1); lnSs_true(f_idx) = lnSs(1); lnSs_true(c_idx) = lnSs(2);

%% Model Inputs and Simulate Data
% Create inputs needed for OHT3D
[inputs] = OHT_create_inputs(well_locs,test_list,domain);
y_fxn = @(params) OHT_run_distribKSs(params, domain, bdrys, inputs, 1);

num_omegas = size(inputs,1);
numobs = size(test_list,1);

%Simulate Synthetic Data - Column vector, where y_oht(1:num_cells) are the
%real phasor component and y_oht(num_cells+1:end) are imaginary phasor
%component
y_oht = y_fxn([lnK_true; lnSs_true]);
A_syn = y_oht(1:numobs);
B_syn = y_oht(numobs+1:2*numobs);

% Calculate head amplitude and phase using numerical steady-periodic head phasors
amp_oht = sqrt(A_syn.^2 + B_syn.^2);
phase_oht = atan2(-B_syn, A_syn);

%% Analytical modeling
soln = 'leaky'; % Specify leaky analytical model for validation

% Constant vertical discretization
T = exp(lnK(1)) * dz * 2;        % Fracture transmissivity (m^2/s)
S = exp(lnSs(1)) * dz * 2;       % Fracture storativity (-)
L = (exp(lnK(2)) / (dz(2)+dz(3))) * 2;  % Confining unit leakance (s^-1)

syn_data(:,1) = (2*pi) ./ test_list(:,1); % Oscillation period (s)
syn_data(:,2) = test_list(:,1);           % Angular Frequency (rad/s)
syn_data(:,3) = test_list(:,3) * 2;       % Max Pumping Rate (m^3/s)
syn_data(:,4) = r(test_list(:,4)-1);      % Inter-well spacing (m)

% Generate steady-periodic head phasors using leaky-confined aquifer analytical model.
for k = 1 : numel(syn_data(:,1))
    y_ana(k,:) = RasSoln(syn_data(k,:), log([T; S; L]), soln);
end

% Calculate head amplitude and phase using analytical head phasors.
amp_ana = sqrt(y_ana(:,1).^2 + y_ana(:,2).^2);
phase_ana = atan2(-y_ana(:,2), y_ana(:,1));

%%Figures
r_list = unique(r);

figure
clf
subplot(1,2,1)
ax = gca;
hold on
for j = 1 : numel(r_list)
    idx = find(syn_data(:,4) == r_list(j));
    plot(P, amp_oht(idx), 'k-', 'LineWidth', 2)
    plot(P, amp_ana(idx), 'k^', 'MarkerSize', 8)
end
xlabel('Period (s)')
ylabel('Amplitude (m)')
ax.XScale = 'log';
ax.FontSize = 30;

subplot(1,2,2)
ax = gca;
hold on
for j = 1 : numel(r_list)
    idx = find(syn_data(:,4) == r_list(j));
    plot(P, phase_oht(idx), 'k-', 'LineWidth', 2)
    plot(P, phase_ana(idx), 'k^', 'MarkerSize', 8)
end
xlabel('Period (s)')
ylabel('Phase (rad)')
ax.XScale = 'log';
ax.FontSize = 30;
set(gcf, 'Position', [100 100 2025 2025/2.6667])

% Uncomment to plot volume slices of model domain
% figure
% clf
% sx = slice(cgrid{1}, cgrid{2}, cgrid{3}, reshape(lnK_true, numy, numx, numz), 0, [], []);
% hold on
% sy = slice(cgrid{1}, cgrid{2}, cgrid{3}, reshape(lnK_true, numy, numx, numz), [], 0, []);
% sz = slice(cgrid{1}, cgrid{2}, cgrid{3}, reshape(lnK_true, numy, numx, numz), [], [], node_list(1));
% sx.FaceColor = 'flat';
% sy.FaceColor = 'flat';
% sz.FaceColor = 'flat';
% axis([-10 10 -10 10 0 0.25])