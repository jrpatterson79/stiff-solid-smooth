% Oscillatory Flow Testing Fracture Characterization
% Numerical Modeling Analysis

% Synthetic modeling experiment to explore fracture apeture heterogeneity as a potential mechanism
% contributing to period dependent fracture flow parameters under
% oscillatory flow conditions. The code uses OHT3D
% (https://github.com/wischydro-cardiff/oscillatory-tomography) to generate
% synthetic data.

% This code generates a mat file that is used to reproduce Figure 2, Figure 6 D-F, and Figure 7 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Stiff, Smooth, and Solid?: Complex Fracture Hydraulic Hydraulics' Imprints on Oscillatory Hydraulic Testing. Submitted to Water Resources Research.

% Code developed by Jeremy Patterson
% Created Dec 2020; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Specify Directory
addpath(genpath('/.../.../')) % Provide location to func_lib subdirectory which contains necessary function files to execute this code

%% Specify Model Run Type
% Choose homogeneous ('homog') or heterogeneous ('heterog') model run
model = 'heterog';
% Specify geostatistical model ('exponential', 'gaussian', 'linear')
geostat_model = 'exponential';

% Set seed for random number generator
num_relz = 1; % Number of realizations for stochastic analysis (set to 100 to reproduce stochastic analysis)
if num_relz == 1
    % Set RNG seed for aperture realizations (Case 1 = 50; Case 2 = 15)
    seed = 15;
    randn('state', seed)
end

%% Load Well Data

% Well locations
well_locs = [0 0;...
             0 -4;...
             -6 0;...
             0 8;...
            10 0];   
num_wells = numel(well_locs(:,1));
%% Forward Model Setup
% Model domain
xmin = -1500; xmax = -xmin; dx = 2;
ymin = xmin; ymax = xmax; dy = dx;

domain = struct('x',[],'y',[],'z',[]);
domain.x = [xmin : dx : xmax];
domain.y = [ymin : dy : ymax];
domain.z = [0 1];

numx = numel(domain.x) - 1;
numy = numel(domain.y) - 1;
num_cells = numx * numy;

% Grid Domain
[coords, cgrid] = plaid_cellcenter_coord(domain);
        
% Coordinates surrounding well field
coords_range = (coords(:,1) > -400 & coords(:,1) < 400 & coords(:,2) > -400 & coords(:,2) < 400);
       
% Specify boundary types and boundary values (x / y constant head
% boundaries, no flux in z)
bdry_types = [1; 1; 1; 1; 0; 0];
bdry_vals = zeros(6,1);
bdry_L = zeros(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals,'leaks',bdry_L);
%% Random field parameters

% Fracture geostatistical properties
ln_aper_mean = log(3e-4);          % Mean aperture (m)
if strcmp(model, 'heterog') == 1
    var_b = 0.3;          % Aperture variance               
elseif strcmp(model, 'homog') == 1
    var_b = 0;            % Aperture variance
end 
L_x = 20; L_y = L_x/4;  % Correlation lengths (m)

% Fluid and Rock Properties
rho = 998.23;           % Fluid density (kg/m^3)
g = 9.81;               % Gravitational acceleration (m/s^2)
mu = 1.0016e-3;         % Fluid dynamic viscosity (Pa s)
RockComp = 1.41e-8;     % Rock compressibility (Pa^-1)
WaterComp = 1 / 2.1e9;  % Fluid compressibility (Pa^-1)
eta_r = 0.15;           % Host rock porosity (-)

% Natural rock fractures contain asperities. The next 5 lines of code take a volume averaging approach to calculate an effective fracture storativity 
% to be applied throughout the fracture. 
Ss_rock = rho * g * (RockComp + (eta_r * WaterComp)); 
Ss_frac = rho * g * WaterComp;

V_frac = Ss_frac * dx(1) * dy(1) * 0.4;
V_rock = Ss_rock * dx(1) * dy(1) * 0.6;
Ss_eff = (V_frac + V_rock) / (dx(1)*dy(1));

% Mean fracture flow parameters
lnT_mean = log((rho*g*exp(ln_aper_mean).^3) ./ (12*mu));
lnS_mean = log(Ss_eff .* exp(ln_aper_mean));
lnD_mean = lnT_mean - lnS_mean;

%% Create Test List
P = logspace(1, 3, 20);  % Stimulation periods [s]
Q_max = 7e-5;            % Peak volumetric flow rate [m^3/s]

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

%% Generate Data / Conduct Inversion
% Test list used for homogeneous parameter estimation
syn_data(:,1) = (2*pi) ./ test_list(:,1); % Oscillation period (s)
syn_data(:,2) = test_list(:,1);           % Angular Frequency (rad/s)
syn_data(:,3) = test_list(:,3) * 2;       % Max Pumping Rate (m^3/s)
syn_data(:,4) = r;                        % Radial Distance (m)

% Generate data covariance matrix
data_err = sqrt(1e-20);        % Expected data measurement error
R_inv = inv(data_err*eye(2));  % Data covariance matrix

% Levenberg-Marquardt inversion parameters
soln = 'confined';
s_init = [lnT_mean; lnS_mean];  % Initial parameter guesses for inversion
delta = [0.1; 0.1];             % Small parameter perturbation to calculate Jacobian
lambda = 1e1;                   % LM initial stabilization parameter
max_iter = 100; 

for relz = 1 : num_relz
    fprintf('Realization %d\n', relz)
    
    % Generate parameter covariance matrix
    distmat_row = dimdist(coords(1,:),coords); % Distance matrix
    max_dist = max(max((distmat_row(:,:,1).^2 + ...
        distmat_row(:,:,2).^2).^.5));
    
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

    % Uncomment these lines to plot aperture realization if desired
    % p = pcolor(cgrid{1}, cgrid{2}, reshape(log10(exp(ln_aper_true)), numy, numx));
    % p.LineStyle = 'none';
    % hold on
    % plot(well_locs(:,1), well_locs(:,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
    % axis([-50 50 -50 50])
    % caxis([-5 -2])
    % pause

    % Convert aperture to flow parameters
    lnT_true = log((rho*g*exp(ln_aper_true).^3) ./ (12*mu));
    lnS_true = log(Ss_eff .* exp(ln_aper_true));
        
    % Generate inputs for OHT3D
    [inputs] = OHT_create_inputs(well_locs,test_list,domain);
    num_omegas = size(inputs,1);
    num_obs = size(test_list,1);
    
    %Simulate Synthetic Data - Column vector, where y_oht(1:num_cells) are the
%real phasor component and y_oht(num_cells+1:end) are imaginary phasor
%component
    y_fxn = @(params) OHT_run_distribKSs(params, domain, bdrys, inputs, 1);
    y_obs = y_fxn([lnT_true; lnS_true]);
    
    % Conduct parameter estimation s_hat = [T S D]
    for i = 1 : num_obs
        y_func = @(s) RasSoln(syn_data(i,:), s, 'confined');
        [pest, ~, out_flag] = Lev_Marq(syn_data(i,:), s_init, [y_obs(i); y_obs(num_obs+i)], R_inv, lambda, delta, soln, max_iter);
        J = jacob(pest, delta, syn_data(i,:), soln);
        
        if out_flag == 0 || s_unc(1) > 1 || s_unc(2) > 1
            [pest, ~, out_flag] = Lev_Marq(syn_data(i,:), pest, [y_obs(i); y_obs(num_obs+i)], R_inv, lambda, delta, soln, max_iter);
            y_opt(i,:) = RasSoln(syn_data(i,:), pest, soln);
            s_hat(i,:) = pest;
        else
            y_opt(i,:) = RasSoln(syn_data(i,:), pest, soln);
            s_hat(i,:) = pest;
        end        
    end
    
    s_opt{relz} = [s_hat s_hat(:,1)-s_hat(:,2)];
    aper(:,relz) = ln_aper_true;
    lnT(:,relz) = lnT_true;
    lnS(:,relz) = lnS_true;
    y_mod{relz} = y_opt;
    synth_data{relz} = syn_data;
end

if num_relz > 1
    % Uncomment and provide a directory location save_dir = /.../.../ to save workspace to a mat file
    % save([save_dir 'stoch_exp_var_aniso.mat'])
elseif num_relz == 1
    % Aperture statistics
    pen_dep = logspace(1,3,2*numel(P));
    aper_geo_mean = zeros(2*numel(P),1);
    aper_ci = zeros(2*numel(P),2);
    for j = 1 : 2*numel(P)
%         pen_dep(j) = sqrt(2*exp(lnD_mean).*P(j));
        idx = sqrt(coords(:,1).^2 + coords(:,2).^2) <= pen_dep(j);
        aper_geo_mean(j) = mean(ln_aper_true(idx));
        aper_ci(j,:) = [prctile(ln_aper_true(idx),5) prctile(ln_aper_true(idx),95)];
    end
    % Uncomment and provide a directory location save_dir = /.../.../ to save workspace to a mat file.
%     save([save_dir 'seed_' num2str(seed) '.mat'])
end