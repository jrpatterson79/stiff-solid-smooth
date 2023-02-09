% Oscillatory Flow Testing Fracture Characterization
% Numerical Modeling Analysis

% Synthetic modeling experiment to explore fracture-host rock fluid exchange as a potential mechanism
% contributing to period dependent fracture flow parameters under
% oscillatory flow conditions. The code uses OHT3D
% (https://github.com/wischydro-cardiff/oscillatory-tomography) to generate
% synthetic data.

% This code generates a mat file that is used to reproduce Figure 6 G-I and Figure 10 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Stiff, Solid, and Smooth?: Complex Fracture Hydraulic Hydraulics Revealed through Oscillatory Flow Interference Testing. Submitted to Water Resources Research.

% Code developed by Jeremy Patterson
% Created Aug 2021; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Specify Directory
file_dir = '/.../.../'; % Add location to COMSOL output csv files ('hydromech_imperm.csv' 'hydromech_perm.csv')
addpath('/.../.../') % Add location to func_lib subdirectory which contains needed function files

%% Gradient Inversion Parameters
% Assumed forward model for parameter estimation
soln = 'confined';

data_err = sqrt(1e-20);         % Assumed data measurement error
R_inv = inv(data_err * eye(2)); % Data error covariance matrix

% Gradient inversion parameters
delta = [0.1; 0.1];  % Jacobian parameter perturbation
lambda = 1e1;        % LM stabilization parameter

% Initial parameter guesses
s_init = [-16 -15;
          -16 -15;
          -16 -16;
          -16 -16;
          -17 -17;
          -17 -18;
          -16 -15;
          -16 -15;
          -16 -16;
          -16 -16;
          -17 -17;
          -17 -18;
          -16 -15;
          -16 -15;
          -16 -16;
          -16 -16;
          -17 -17;
          -17 -17;
          -16 -14;
          -16 -15;
          -16 -15;
          -16 -16;
          -17 -17;
          -17 -17;
          -16 -14;
          -16 -15;
          -16 -15;
          -16 -15;
          -16 -16;
          -17 -17;
          -16 -14;
          -16 -15;
          -16 -15;
          -16 -15;
          -16 -16;
          -17 -17;
          -16 -14;
          -16 -14;
          -16 -15;
          -16 -15;
          -16 -15;
          -16 -16;
          -16 -13;
          -16 -14;
          -16 -15;
          -16 -15;
          -16 -15;
          -16 -16;
          -16 -13;
          -16 -14;
          -16 -14;
          -16 -15;
          -16 -15;
          -16 -15;
          -16 -13;
          -16 -14;
          -16 -14;
          -16 -15;
          -16 -15;
          -16 -15;
          -16 -13;
          -16 -14;
          -16 -14;
          -16 -14;
          -16 -15;
          -16 -15;
          -16 -13;
          -16 -13;
          -16 -14;
          -16 -14;
          -16 -14;
          -16 -15;
          -16 -12;
          -16 -13;
          -16 -14;
          -16 -14;
          -16 -14;
          -16 -15;
          -16 -12;
          -16 -13;
          -16 -13;
          -16 -14;
          -16 -14;
          -16 -14;
          -16 -12;
          -16 -13;
          -16 -13;
          -16 -14;
          -16 -14;
          -16 -14;
          -16 -12;
          -16 -13;
          -16 -13;
          -16 -13;
          -16 -14;
          -16 -14;
          -17 -12;
          -16 -12;
          -16 -13;
          -16 -13;
          -16 -13;
          -16 -14;
          -16 -11;
          -16 -12;
          -16 -13;
          -16 -13;
          -16 -13;
          -16 -14;
          -16 -11;
          -16 -12;
          -16 -12;
          -16 -13;
          -16 -13;
          -16 -13;
          -16 -11;
          -16 -12;
          -16 -12;
          -16 -13;
          -16 -13;
          -16 -13;
          -17 -11;
          -16 -12;
          -16 -12;
          -16 -12;
          -16 -13;
          -16 -13];

%% Impermeable Host Rock Analysis
% Load data
file = 'hydromech_imperm.csv';
data = readtable([file_dir file]);

% Stimulation Periods (s)
P = unique(data.P_s_);
num_p = numel(P);

% Well Distances
r = [4; 6; 8; 10; 12; 16];
num_r = numel(r);

Q_max = 7e-5; % Max pumping rate (m^3/s)
% Generate Test List
test_list = [];
for i = 1 : num_p
    for j = 1 : num_r
        test_list = [test_list;
                     P(i) (2*pi)/P(i) Q_max r(j) j];
    end
end
num_tests = numel(test_list(:,1));

% Observed Head (m)
conv_fact = 1.019977e-4; % Pa to m head conversion factor
h_obs = [data.Obs1_Hyd data.Obs2_Hyd data.Obs3_Hyd data.Obs4_Hyd data.Obs5_Hyd data.Obs6_Hyd] * conv_fact;
numobs = numel(h_obs(1,:));

% Observed fracture displacement (m)
strain_obs = [data.Obs1_Mech data.Obs2_Mech data.Obs3_Mech data.Obs4_Mech data.Obs5_Mech data.Obs6_Mech];
num_strain = numel(strain_obs);

% Extract pressure and strain phasors; Homogeneous parameter estimation
press_phasor_imperm = zeros(num_tests,1); strain_phasor_imperm = zeros(num_tests,1);
for i = 1 : num_tests
    % Find steady-periodic portion of signal
    idx = data.time_s_ >= 5*test_list(i,1) & data.P_s_ == test_list(i,1);
    
    % Extract steady-periodic head and displacement phasors
    [~, press_phasor_imperm(i)] = periodic_LS_fit(data.time_s_(idx), h_obs(idx,test_list(i,5)), test_list(i,1));
    [~, strain_phasor_imperm(i)] = periodic_LS_fit(data.time_s_(idx), strain_obs(idx, test_list(i,5)), test_list(i,1));
    
    % Reconstructed pressure and displacement signals
    press_recon = (real(press_phasor_imperm(i)) .* cos(test_list(i,2) .* data.time_s_(idx))) -...
                  (imag(press_phasor_imperm(i)) .* sin(test_list(i,2) .* data.time_s_(idx)));
    strain_recon = (real(strain_phasor_imperm(i)) .* cos(test_list(i,2) .* data.time_s_(idx))) -...
                   (imag(strain_phasor_imperm(i)) .* sin(test_list(i,2) .* data.time_s_(idx)));
    
    % Conduct parameter estimation
    [s_hat_imperm(i,:), ~, flag_imperm(i)] = Lev_Marq(test_list(i,:), s_init(i,:)', [real(press_phasor_imperm(i)); imag(press_phasor_imperm(i))], R_inv, lambda, delta, soln);
    y_opt_imperm(i,:) = RasSoln(test_list(i,:), s_hat_imperm(i,:), soln);
end
s_hat_imperm(:,3) = s_hat_imperm(:,1) - s_hat_imperm(:,2);
rmse_imperm = sqrt(mean([real(press_phasor_imperm); imag(press_phasor_imperm)] - [y_opt_imperm(:,1); y_opt_imperm(:,2)]).^2);

%% Permeable Host Rock Analysis (k = 3e-15)
% Load data
file = 'hydromech_perm.csv';
data_perm = readtable([file_dir file]);

% Observed Head (m)
h_perm = [data_perm.Obs1_Hyd data_perm.Obs2_Hyd data_perm.Obs3_Hyd data_perm.Obs4_Hyd data_perm.Obs5_Hyd data_perm.Obs6_Hyd] * conv_fact;

% Observed Deformation (m)
strain_perm = [data_perm.Obs1_Mech data_perm.Obs2_Mech data_perm.Obs3_Mech data_perm.Obs4_Mech data_perm.Obs5_Mech data_perm.Obs6_Mech];

% Extract Fourier Coefficients
press_phasor_perm = zeros(num_tests,1); strain_phasor_perm = zeros(num_tests,1);
for i = 1 : num_tests
    % Find steady-periodic portion of signal
    idx = data_perm.time_s_ >= 5*test_list(i,1) & data_perm.P_s_ == test_list(i,1);
    
    % Extract steady-periodic pressure and fracture displacement phasors
    [~, press_phasor_perm(i)] = periodic_LS_fit(data_perm.time_s_(idx), h_perm(idx,test_list(i,5)), test_list(i,1));
    [~, strain_phasor_perm(i)] = periodic_LS_fit(data_perm.time_s_(idx), strain_perm(idx, test_list(i,5)), test_list(i,1));
    
    % Reconstructed pressure and displacement signals
    press_recon_perm = (real(press_phasor_perm(i)) .* cos(test_list(i,2) .* data_perm.time_s_(idx))) -...
                       (imag(press_phasor_perm(i)) .* sin(test_list(i,2) .* data_perm.time_s_(idx)));
    strain_recon_perm = (real(strain_phasor_perm(i)) .* cos(test_list(i,2) .* data_perm.time_s_(idx))) -...
                        (imag(strain_phasor_perm(i)) .* sin(test_list(i,2) .* data_perm.time_s_(idx)));
    
    % Conduct parameter estimation
    [s_hat_perm(i,:), ~, flag_perm(i)] = Lev_Marq(test_list(i,:), s_init(i,:)', [real(press_phasor_perm(i)); imag(press_phasor_perm(i))], R_inv, lambda, delta, soln);
    y_opt_perm(i,:) = RasSoln(test_list(i,:), s_hat_imperm(i,:), soln);
end
s_hat_perm(:,3) = s_hat_perm(:,1) - s_hat_perm(:,2);

%% Permeable Host Rock Analysis (k = 3e-13)
% Load data
file = 'hydromech_k13.csv';
data_13 = readtable([file_dir file]);

% Observed Head (m)
h_13 = [data_13.Obs1_Hyd data_13.Obs2_Hyd data_13.Obs3_Hyd data_13.Obs4_Hyd data_13.Obs5_Hyd data_13.Obs6_Hyd] * conv_fact;

% Observed Deformation (m)
strain_13 = [data_13.Obs1_Mech data_13.Obs2_Mech data_13.Obs3_Mech data_13.Obs4_Mech data_13.Obs5_Mech data_13.Obs6_Mech];

% Extract phasor Coefficients
press_phasor_13 = zeros(num_tests,1); strain_phasor_13 = zeros(num_tests,1);
for i = 1 : num_tests
    % Find steady-periodic portion of signal
    idx = data_13.time_s_ >= 5*test_list(i,1) & data_13.P_s_ == test_list(i,1);
    
    % Extract steady-periodic head and fracture displacement phasors
    [~, press_phasor_13(i)] = periodic_LS_fit(data_13.time_s_(idx), h_13(idx,test_list(i,5)), test_list(i,1));
    [~, strain_phasor_13(i)] = periodic_LS_fit(data_13.time_s_(idx), strain_13(idx, test_list(i,5)), test_list(i,1));
    
    % Reconstructed pressure and displacement signals
    press_recon_13 = (real(press_phasor_13(i)) .* cos(test_list(i,2) .* data_13.time_s_(idx))) -...
                       (imag(press_phasor_13(i)) .* sin(test_list(i,2) .* data_13.time_s_(idx)));
    strain_recon_13 = (real(strain_phasor_13(i)) .* cos(test_list(i,2) .* data_13.time_s_(idx))) -...
                        (imag(strain_phasor_13(i)) .* sin(test_list(i,2) .* data_13.time_s_(idx)));
    
    % Conduct parameter estimation
    [s_hat_13(i,:), ~, flag_13(i)] = Lev_Marq(test_list(i,:), s_init(i,:)', [real(press_phasor_13(i)); imag(press_phasor_13(i))], R_inv, lambda, delta, soln);
    y_opt_13(i,:) = RasSoln(test_list(i,:), s_hat_13(i,:), soln);
end
s_hat_13(:,3) = s_hat_13(:,1) - s_hat_13(:,2);

% Uncomment to save workspace to mat file if desired
% save_name = 'hydromechanical_pest_results.mat';
% save([file_dir save_name])