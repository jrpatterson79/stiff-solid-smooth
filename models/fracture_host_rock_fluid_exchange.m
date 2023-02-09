% Oscillatory Flow Testing Fracture Characterization
% Numerical Modeling Analysis

% Synthetic modeling experiment to explore fracture-host rock fluid exchange as a potential mechanism
% contributing to period dependent fracture flow parameters under
% oscillatory flow conditions. The code uses OHT3D
% (https://github.com/wischydro-cardiff/oscillatory-tomography) to generate
% synthetic data.

% This code generates a mat file that is used to reproduce Figure 6 G-I, and Figure 8 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Stiff, Solid, and Smooth?: Complex Fracture Hydraulic Hydraulics Revealed through Oscillatory Flow Interference Testing. Submitted to Water Resources Research.

% Code developed by Jeremy Patterson
% Created Dec 2020; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Specify Directory
save_dir = '/Users/jpatt/My Drive/PhD/WRR_FreqDepend/MAT_Files/Fluid_Exchange/';
fig_dir = '/Users/jpatt/My Drive/PhD/Drafts/Patterson_Cardiff_WRR_2023/Figures/';
addpath(genpath('/Users/jpatt/My Drive/PhD/WRR_FreqDepend/Code_Sandbox/'))

%% Specify Model Run Type
% Choose homogeneous ('homog') or heterogeneous ('heterog') model run
model = 'homog';   

% Specify geostatistical model ('exponential', 'gaussian', 'linear')
geostat_model = 'exponential';

% Set RNG seed for aperture realizations (Case 1 = 50; Case 2 = 15)
seed = 50;
randn('state', seed)

% Specify whether to use confined or leaky analytical solution
soln = 'confined';
%% Forward Model Setup

% Specify Domain Geometry
xmin = -400; xmax = -xmin; dx = 2;
ymin = xmin; ymax = xmax; dy = dx;
dz = 5e-2; zmin = 0; zmax = 5*dz;

domain = struct('x',[],'y',[],'z',[]);
domain.x = [xmin : dx : xmax];
domain.y = [ymin : dy : ymax];
domain.z = [zmin: dz : zmax];

numx = numel(domain.x) - 1;
numy = numel(domain.y) - 1;
numz = numel(domain.z) - 1;
num_cells = numx * numy * numz;

% Specify well locations
elev = dz / 2; 
well_locs = [0 0 elev;...
             0 -4 elev;...
             -6 0 elev;...
             0 8 elev;...
            10 0 elev];
num_wells = numel(well_locs(:,1));

% Grid Domain
[coords, cgrid] = plaid_cellcenter_coord(domain);
node_list = unique(coords(:,3));
f_idx = find(coords(:,3) == node_list(1));     % Locate fracture grid cells in the model domain
fx_coords = [coords(f_idx,1) coords(f_idx,2)]; 

% Specify boundary types and boundary values (x / y constant head
% boundaries, no flux in z)
bdry_types = [2; 2; 2; 2; 0; 0];
bdry_vals = zeros(6,1);
bdry_L = 1e-5*ones(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals,'leaks',bdry_L);

%% Create Test List
P = logspace(1, 3, 20); % Oscilltion period [s]
Q_max = 7e-5;           % Max pumping rate [m^3/s]

%Test list for OHT input [omega pump_well_index peak_flow_rate obs_well_index]
test_list = [];
r = [];
for i = 1:numel(P)
    for j = 1 : num_wells
        for k = j+1:num_wells
            r  = [r; 
                  sqrt((well_locs(j,1)-well_locs(k,1)).^2 + (well_locs(j,2)-well_locs(k,2)).^2)];
            test_list = [test_list; ...
                         (2*pi)/P(i) j Q_max k];
        end
    end
end
num_tests = numel(test_list(:,1));

%% Specify Model Parameters
% Fracture geostatistical properties
ln_aper_mean = log(3e-4); % Mean aperture (m)

if strcmp(model, 'heterog') == 1
    var_b = 0.3;          % Aperture variance
elseif strcmp(model, 'homog') == 1
    var_b = 0;            % Aperture variance
end
% Correlation lengths (m)
L_x = 20; L_y = L_x / 4;                 

% Distance matrix for geostatistical realization
distmat_row = dimdist(fx_coords(1,:), fx_coords);
max_dist = max(max((distmat_row(:,:,1).^2 + ...
    distmat_row(:,:,2).^2).^.5));

% Generate the correlation matrix and geostatistical realizations based on
% the chosen variogram model
if strcmp(geostat_model, 'exponential') == 1
    corr_row = exp(-(...
        (distmat_row(:,:,1)./L_x).^2 + ...
        (distmat_row(:,:,2)./L_y).^2).^.5);
    [corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[numy numx]);
    ln_aper_true = corr_relz(:,1).*var_b.^.5 + ln_aper_mean;
    
elseif strcmp(geostat_model, 'gaussian') == 1
    corr_row = exp(-(...
        (distmat_row(:,:,1)./L_x).^2 + ...
        (distmat_row(:,:,2)./L_y).^2));
    [corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[numy numx]);
    ln_aper_true = corr_relz(:,1).*var_b.^.5 + ln_aper_mean;
    
elseif strcmp(geostat_model, 'linear') == 1
    corr_row = (1/max_dist) .*...
        (max_dist - ((distmat_row(:,:,1).^2 + ...
        distmat_row(:,:,2).^2).^.5));
    [corr_relz] = toepmat_vector_math(corr_row,'r',[],2,[numy numx]);
    ln_aper_true = corr_relz(:,1).*(1.8622*var_b).^.5 + ln_aper_mean;
    
else
    error('Pick a valid variogram model')
end

% Fluid and Rock Properties
rho = 998.23;           % Fluid density (kg/m^3)
g = 9.81;               % Gravitational acceleration (m/s^2)
mu = 1.0016e-3;         % Fluid dynamic viscosity (Pa s)
RockComp = 1.41e-8;     % Rock compressibility (Pa^-1)
WaterComp = 1 / 2.1e9;  % Fluid compressibility (Pa^-1)
eta_r = 0.15;           % Host rock porosity (-)

% Fracture properties
T_frac = ((rho*g)/(12 * mu)) .* exp(ln_aper_true).^3; 
K_frac = log(T_frac ./ (dz*2)); % Fracture hydraulic conductivity
K_rock = log(3e-8);             % Host rock hydraulic conductivity

Ss_rock = rho * g * (RockComp + (eta * WaterComp)); % Host rock specific storage
% Natural rock fractures contain asperities. The next 5 lines of code take a volume averaging approach to calculate an effective fracture storativity 
% to be applied throughout the fracture. 
Ss_f = rho * g * WaterComp;
V_frac = Ss_f * dx * dy * 0.4;
V_rock = Ss_rock * dx * dy * 0.6;
Ss_eff = (V_frac + V_rock) / (dx*dy);
S_frac = Ss_eff .* exp(ln_aper_true);
Ss_frac = log(S_frac ./ (dz*2)); 

lnK_true = K_rock * ones(num_cells,1);   lnK_true(f_idx) = K_frac;
lnSs_true = log(Ss_rock) * ones(num_cells,1); lnSs_true(f_idx) = Ss_frac;

%% Model Inputs and Simulate Data
% Create inputs needed for OHT3D
[inputs] = OHT_create_inputs(well_locs,test_list,domain);
y_fxn = @(params) OHT_run_distribKSs(params, domain, bdrys, inputs, 1);

num_omegas = size(inputs,1);
num_obs = size(test_list,1);

%Simulate Synthetic Data - Column vector, where y_oht(1:num_cells) are the
%real phasor component and y_oht(num_cells+1:end) are imaginary phasor
%component
y_oht = y_fxn([lnK_true; lnSs_true]);

%% Homogeneous Parameter Estimation
% Test list used for homogeneous parameter estimation
syn_data(:,1) = (2*pi) ./ test_list(:,1); % Oscillation period (s)
syn_data(:,2) = test_list(:,1);           % Angular Frequency (rad/s)
syn_data(:,3) = test_list(:,3) * 2;       % Max Pumping Rate (m^3/s)
syn_data(:,4) = r;                        % Radial Distance (m)

data_err = sqrt(1e-20);          % Assumed data measurement error [m] 
R_inv = inv(data_err * eye(2));  % Data error covariance matrix

% Initial Inversion Parameters
lambda_init = 1e1;  % LM initial stabilization parameter
delta = [0.1; 0.1]; % Parameter perturbation parameter for Jacobian
s_init = [-16 -15;
          -16 -15;
          -15 -15;
          -14 -15;
          -15 -15;
          -14 -14;
          -14 -15;
          -14 -15;
          -13 -14;
          -13 -14;
          -17 -15;
          -16 -15;
          -15 -15;
          -14 -14;
          -15 -15;
          -14 -14;
          -14 -14;
          -14 -14;
          -12 -13;
          -13 -14;
          -16 -14;
          -16 -15;
          -15 -15;
          -15 -15;
          -15 -15;
          -14 -14;
          -14 -14;
          -14 -14;
          -13 -14;
          -14 -14;
          -16 -14;
          -16 -15;
          -15 -14;
          -15 -15;
          -16 -15;
          -14 -14;
          -15 -14;
          -15 -14;
          -13 -14;
          -14 -14;
          -17 -14;
          -16 -14;
          -15 -14;
          -15 -14;
          -16 -14;
          -14 -14;
          -15 -14;
          -15 -14;
          -13 -13;
          -14 -14;
          -17 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -16 -14;
          -15 -14;
          -15 -14;
          -15 -14;
          -14 -14;
          -15 -14;
          -17 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -16 -14;
          -15 -14;
          -15 -14;
          -15 -14;
          -14 -14;
          -15 -14;
          -17 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -15 -14;
          -14 -14;
          -15 -14;
          -17 -14;
          -16 -14;
          -16 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -15 -14;
          -16 -14;
          -15 -14;
          -15 -14;
          -17 -14;
          -17 -14;
          -16 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -15 -14;
          -17 -13;
          -17 -14;
          -16 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -15 -14;
          -17 -13;
          -17 -13;
          -16 -14;
          -16 -14;
          -16 -14;
          -16 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -16 -14;
          -17 -13;
          -17 -13;
          -16 -14;
          -16 -14;
          -17 -14;
          -16 -14;
          -16 -14;
          -16 -14;
          -15 -14;
          -16 -14;
          -17 -13;
          -17 -13;
          -16 -13;
          -16 -14;
          -16 -13;
          -16 -14;
          -16 -13;
          -16 -14;
          -15 -14;
          -16 -14;
          -17 -13;
          -17 -13;
          -16 -13;
          -16 -14;
          -17 -14;
          -16 -14;
          -16 -14;
          -16 -13;
          -16 -14;
          -16 -14;
          -17 -12;
          -17 -13;
          -17 -13;
          -16 -13;
          -17 -13;
          -16 -13;
          -16 -13;
          -16 -13;
          -16 -14;
          -16 -13;
          -17 -12;
          -17 -13;
          -17 -13;
          -17 -14;
          -17 -13;
          -17 -14;
          -17 -14;
          -17 -13;
          -16 -14;
          -16 -13;
          -17 -12;
          -17 -13;
          -17 -13;
          -17 -13;
          -17 -13;
          -17 -14;
          -17 -13;
          -17 -13;
          -16 -13;
          -17 -14;
          -18 -12;
          -17 -12;
          -17 -13;
          -17 -13;
          -17 -13;
          -17 -13;
          -17 -13;
          -17 -13;
          -16 -14;
          -17 -14;
          -18 -12;
          -17 -12;
          -17 -13;
          -17 -13;
          -17 -13;
          -17 -13;
          -17 -13;
          -17 -13;
          -17 -13;
          -17 -13];

% Uncomment for realization 2 (seed = 15)
% s_init = [-16 -14;
%           -15 -14;
%           -15 -15;
%           -14 -14;
%           -13 -13;
%           -11 -12;
%           -11 -12;
%           -14 -14;
%           -12 -13;
%           -14 -15;
%           -16 -14;
%           -15 -14;
%           -15 -15;
%           -14 -14;
%           -13 -12;
%           -12 -12;
%           -11 -11;
%           -14 -14;
%           -12 -13;
%           -14 -14;
%           -16 -14;
%           -15 -14;
%           -15 -14;
%           -14 -14;
%           -14 -13;
%           -12 -12;
%           -11 -11;
%           -14 -14;
%           -12 -13;
%           -14 -14;
%           -16 -14;
%           -16 -14;
%           -15 -14;
%           -14 -14;
%           -14 -13;
%           -12 -12;
%           -11 -11;
%           -15 -14;
%           -13 -13;
%           -14 -14;
%           -16 -14;
%           -16 -14;
%           -15 -14;
%           -15 -14;
%           -14 -13;
%           -13 -12;
%           -12 -11;
%           -15 -14;
%           -13 -13;
%           -15 -14;
%           -17 -14;
%           -16 -14;
%           -16 -14;
%           -15 -14;
%           -15 -13;
%           -13 -12;
%           -12 -11;
%           -15 -14;
%           -13 -13;
%           -15 -14;
%           -17 -14;
%           -16 -14;
%           -16 -14;
%           -15 -14;
%           -15 -13;
%           -13 -12;
%           -13 -11;
%           -15 -14;
%           -14 -13;
%           -15 -14;
%           -17 -14;
%           -16 -14;
%           -16 -14;
%           -15 -14;
%           -15 -13;
%           -14 -12;
%           -13 -11;
%           -16 -14;
%           -14 -13;
%           -15 -14;
%           -17 -13;
%           -17 -14;
%           -16 -14;
%           -16 -14;
%           -16 -13;
%           -14 -12;
%           -13 -11;
%           -16 -14;
%           -14 -13;
%           -15 -14;
%           -17 -13;
%           -17 -14;
%           -16 -14;
%           -16 -14;
%           -16 -13;
%           -14 -12;
%           -14 -12;
%           -16 -14;
%           -15 -14;
%           -16 -14;
%           -17 -13;
%           -17 -14;
%           -17 -14;
%           -16 -14;
%           -16 -13;
%           -15 -13;
%           -14 -12;
%           -16 -14;
%           -15 -14;
%           -16 -14;
%           -17 -13;
%           -17 -13;
%           -17 -14;
%           -16 -14;
%           -16 -13;
%           -15 -13;
%           -14 -12;
%           -16 -14;
%           -15 -13;
%           -16 -14;
%           -18 -13;
%           -17 -13;
%           -17 -14;
%           -16 -13;
%           -17 -13;
%           -15 -13;
%           -15 -12;
%           -17 -14;
%           -15 -13;
%           -16 -14;
%           -18 -13;
%           -17 -13;
%           -17 -14;
%           -16 -13;
%           -17 -13;
%           -15 -12;
%           -15 -12;
%           -17 -14;
%           -15 -13;
%           -16 -14;
%           -18 -13;
%           -17 -13;
%           -17 -14;
%           -16 -13;
%           -17 -13;
%           -15 -12;
%           -15 -12;
%           -17 -14;
%           -16 -13;
%           -16 -13;
%           -18 -13;
%           -17 -13;
%           -17 -13;
%           -17 -13;
%           -17 -13;
%           -16 -13;
%           -16 -12;
%           -17 -14;
%           -16 -13;
%           -17 -14;
%           -18 -13;
%           -18 -13;
%           -17 -13;
%           -17 -13;
%           -18 -13;
%           -16 -12;
%           -16 -12;
%           -17 -13;
%           -16 -13;
%           -17 -14;
%           -19 -13;
%           -18 -13;
%           -17 -13;
%           -17 -13;
%           -18 -13;
%           -17 -13;
%           -16 -12;
%           -17 -13;
%           -16 -13;
%           -17 -14;
%           -19 -13;
%           -18 -13;
%           -18 -13;
%           -17 -13;
%           -18 -13;
%           -17 -13;
%           -17 -13;
%           -17 -13;
%           -17 -14;
%           -17 -13;
%           -19 -13;
%           -18 -13;
%           -18 -13;
%           -17 -13;
%           -18 -13;
%           -17 -13;
%           -17 -12;
%           -18 -13;
%           -17 -13;
%           -17 -13];

% Conduct parameter estimation
for i = 1 : num_obs
    y_func = @(s) RasSoln(syn_data(i,:), s, 'confined');
    [s_opt(i,:), s_step, flag] = Lev_Marq(syn_data(i,:), s_init(i,:)', [y_oht(i); y_oht(num_obs+i)], R_inv, lambda_init, delta, soln);
     y_opt(i,:) = y_func(s_opt(i,:));    
end
s_opt(:,3) = s_opt(:,1) - s_opt(:,2);

% Uncomment to save workspace as a mat file
% if strcmp(model, 'heterog') == 1
%     save([save_dir 'seed_' num2str(seed) '_3D.mat'])
% elseif strcmp(model, 'homog') == 1
%     save([save_dir 'homog_3D.mat'])
% end

%% Figures
% Uncomment to plot volume slices of model domain
% figure
% clf
% ax = gca;
% sx = slice(cgrid{1}, cgrid{2}, cgrid{3}, reshape(lnK_true, numy, numx, numz), 15, [], []);
% hold on
% sy = slice(cgrid{1}, cgrid{2}, cgrid{3}, reshape(lnK_true, numy, numx, numz), [], 15, []);
% sz = slice(cgrid{1}, cgrid{2}, cgrid{3}, reshape(lnK_true, numy, numx, numz), [], [], node_list(1));
% sx.EdgeColor = 'none'; sx.FaceAlpha = 0.7;
% sy.EdgeColor = 'none'; sy.FaceAlpha = 0.7;
% sz.EdgeColor = 'none';
% sx.FaceColor = 'flat';
% sy.FaceColor = 'flat';
% sz.FaceColor = 'flat';
% plot3(well_locs(1,1), well_locs(1,2), well_locs(1,3), 'ko',...
%       'MarkerFaceColor', 'w', 'MarkerSize', 12)
% plot3(well_locs(2:end,1), well_locs(2:end,2), well_locs(2:end,3), 'ko',...
%       'MarkerFaceColor', 'r', 'MarkerSize', 12)
% axis([-50 50 -50 50 node_list(1) zmax])
% ax.XTick = [-50:25:50];
% ax.YTick = [-50:25:50];
% ax.ZTick = [0:0.1:0.2];
% view(-43, 75)
% c1 = colorbar;
% caxis([-15 -3])
% c1.Ticks = [-15:3:-3];
% c1.Label.String = 'ln(K [m/s])';
% c1.FontSize = 24;
% xlabel('X (m)')
% ylabel('Y (m)')
% zlabel('Z (m)')
% ax.FontSize = 28;
% set(gcf, 'Position', [100 100 975 975/1.5])
