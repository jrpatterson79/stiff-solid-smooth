% Oscillatory Flow Testing Fracture Characterization
% Numerical Modeling Analysis

% This code validates the 2D numercial model used to explore the impacts of
% aperture heterogeneity on oscillatory flow signals. The numerical model is
% validated using the fully confined analytical model developed by
% Rasmussen et al. (2003).

% Code developed by Jeremy Patterson
% Created Dec 2018; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Specify Directory
addpath(genpath('/.../.../')) % Provide location to func_lib subdirectory that contains the needed function files

%% Forward Model
%Specify domain 
xmin = -1500; xmax = -xmin; dx = [1 2];
ymin = xmin; ymax = xmax; dy = dx;

domain = struct('x',[],'y',[],'z',[]);
domain.x = [xmin : dx(1) : xmax];
domain.y = [ymin : dy(1) : ymax];
domain.z = [0 1];

domain2 = struct('x',[],'y',[],'z',[]);
domain2.x = [xmin : dx(2) : xmax];
domain2.y = [ymin : dy(2) : ymax];
domain2.z = [0 1];

numx = numel(domain.x) - 1;
numy = numel(domain.y) - 1;
num_cells = numx * numy;

numx_2 = numel(domain2.x) - 1;
numy_2 = numel(domain2.y) - 1;
num_cells2 = numx_2 * numy_2;

well_locs = [0 0;...
             0 -4;...
             -6 0;...
             0 8;...
            10 0];
num_wells = numel(well_locs(:,1));

% Grid Domain
[coords, cgrid] = plaid_cellcenter_coord(domain);
[c2, cg2] = plaid_cellcenter_coord(domain2);

% Specify boundary types and boundary values (x / y constant head
% boundaries, no flux in z)
bdry_types = [1; 1; 1; 1; 0; 0];
bdry_vals = zeros(6,1);
bdry_L = zeros(6,1);
bdrys = struct('types',bdry_types,'vals',bdry_vals,'leaks',bdry_L);

%% Create Test List
P = logspace(1, 3, 5);
Q_max = 7e-5;

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

%% 1 m discretization validation.
% Fluid and Rock Properties
rho = 998.23;           % Fluid density (kg/m^3)
g = 9.81;               % Gravitational acceleration (m/s^2)
mu = 1.0016e-3;         % Fluid dynamic viscosity (Pa s)
RockComp = 1.41e-8;     % Rock compressibility (Pa^-1)
WaterComp = 1 / 2.1e9;  % Fluid compressibility (Pa^-1)
eta_r = 0.15;           % Host rock porosity (-)

ln_aper = log(3e-4);                           % Mean fracture aperture (m)
lnT = log((rho*g)/(12*mu) * exp(ln_aper).^3);  % Fracture transmissivity (m^2/s)

% Natural rock fractures contain asperities. The next 5 lines of code take a volume averaging approach to calculate an effective fracture storativity 
% to be applied throughout the fracture. 
Ss_rock = rho * g * (RockComp + (eta * WaterComp)); 
Ss_frac = rho * g * WaterComp;
V_frac = Ss_frac * dx(1) * dy(1) * 0.4;
V_rock = Ss_rock * dx(1) * dy(1) * 0.6;
Ss_eff = (V_frac + V_rock) / (dx(1)*dy(1));
lnS = log(Ss_eff * exp(ln_aper));            % Fracture storativity (-)

tic
% Generate the inputs needed for OHT3D
[inputs] = OHT_create_inputs(well_locs,test_list,domain);
input_time = toc

y_fxn = @(params) OHT_run_distribKSs(params, domain, bdrys, inputs, 1);

num_omegas = size(inputs,1);
num_obs = size(test_list,1);

tic
%Simulate Synthetic Data - Column vector, where y_oht(1:num_cells) are the
%real phasor component and y_oht(num_cells+1:end) are imaginary phasor
%component
y_oht = y_fxn([lnT*ones(num_cells,1); lnS*ones(num_cells,1)]);
phasor_time = toc

A_syn = y_oht(1:num_obs);
B_syn = y_oht(num_obs+1:2*num_obs);

% Use head phasors to calculate head amplitude and phase
amp_oht = sqrt(A_syn.^2 + B_syn.^2);
phase_oht = atan2(-B_syn, A_syn);

%% Analytical modeling for validation
soln = 'confined';

syn_data(:,1) = (2*pi) ./ test_list(:,1); %Stimulation period (s)
syn_data(:,2) = test_list(:,1);           %Angular Frequency (rad/s)
syn_data(:,3) = test_list(:,3);           %Max Pumping Rate (m^3/s)
syn_data(:,4) = r;                        %Radial Distance (m)

% Calculate head phasor using Rasmussen et al. (2003) analytical solution
for k = 1 : numel(syn_data(:,1))
    y_ana(k,:) = RasSoln(syn_data(k,:), [lnT; lnS], soln);
end

amp_ana = sqrt(y_ana(:,1).^2 + y_ana(:,2).^2);
phase_ana = atan2(-y_ana(:,2), y_ana(:,1));

%% 2 m discretization validation
V_frac2 = Ss_frac * dx(2) * dy(2) * 0.4;
V_rock2 = Ss_rock * dx(2) * dy(2) * 0.6;
Ss_eff2 = (V_frac2 + V_rock2) / (dx(2)*dy(2));
lnS2 = log(Ss_eff2 * exp(ln_aper));

tic
[inputs2] = OHT_create_inputs(well_locs,test_list,domain2);
input_time = toc

y_fxn2 = @(params) OHT_run_distribKSs(params, domain2, bdrys, inputs2, 1);

tic
y_oht2 = y_fxn2([lnT*ones(num_cells2,1); lnS2*ones(num_cells2,1)]);
phasor_time = toc
A_syn2 = y_oht2(1:num_obs);
B_syn2 = y_oht2(num_obs+1:2*num_obs);

amp_oht2 = sqrt(A_syn2.^2 + B_syn2.^2);
phase_oht2 = atan2(-B_syn2, A_syn2);

for k = 1 : numel(syn_data(:,1))
    y_ana2(k,:) = RasSoln(syn_data(k,:), [lnT; lnS2], soln);
end

amp_ana2 = sqrt(y_ana2(:,1).^2 + y_ana2(:,2).^2);
phase_ana2 = atan2(-y_ana2(:,2), y_ana2(:,1));

% Model misfit
misfit = mean((y_oht-y_oht2).^2);

%% Figures
% Uncomment to plot aperture field in the model domain
% figure
% clf
% ax = gca;
% p = pcolor(cgrid{1}, cgrid{2}, ln_aper*ones(numy,numx));
% p.LineStyle = 'none';
% hold on
% plot(well_locs(1,1), well_locs(1,2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 12)
% plot(well_locs(2:end,1), well_locs(2:end,2), 'ko', 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'MarkerSize', 12)
% text(well_locs(2:end,1), well_locs(2:end,2), {'O1', 'O2', 'O3', 'O4'}, 'FontSize', 28, 'FontWeight', 'bold',...
%      'VerticalAlignment', 'bottom', 'HorizontalAlignment','center')
% axis([-15 15 -15 15])
% xlabel('X (m)')
% ylabel('Y (m)')
% ax.FontSize = 32;
% c = colorbar;
% caxis([-9 -7])
% c.Label.String = 'ln(Aperture [m])';
% c.FontSize = 28;
% set(gcf, 'Position', [100 100 900 900/1.1])

r_list = unique(r);
figure 
clf
subplot(1,2,1)
ax = gca;
hold on
for j = 1 : numel(r_list)
    idx = find(syn_data(:,4) == r_list(j));
    plot(P, amp_ana(idx), 'k^', 'MarkerSize', 12, 'LineWidth', 2)
    plot(P, amp_oht(idx), 'k-', 'LineWidth', 2)
    plot(P, amp_oht2(idx), 'r--', 'LineWidth', 2)
end
xlabel('Period (s)')
ylabel('Amplitude (m)')
ax.XScale = 'log';
ax.FontSize = 30;
l = legend('Analytical', 'Numerical - 1 m grid', 'Numerical - 2 m grid');
l.Location = 'northwest';
l.FontSize = 30;
text(775, 2.45, 'A', 'FontSize', 30, 'FontWeight', 'Bold')

subplot(1,2,2)
ax = gca;
hold on
for j = 1 : numel(r_list)
    idx = find(syn_data(:,4) == r_list(j));
    plot(P, phase_ana(idx), 'k^', 'MarkerSize', 12, 'LineWidth', 2)
    plot(P, phase_oht(idx), 'k-', 'LineWidth', 2)
    plot(P, phase_oht2(idx), 'r--', 'LineWidth', 2)
end
xlabel('Period (s)')
ylabel('Phase (rad)')
ax.XScale = 'log';
ax.FontSize = 30;
text(800, 0.49, 'B', 'FontSize', 30, 'FontWeight', 'Bold')
set(gcf, 'Position', [100 100 1900 600])

figure
clf
ax = gca;
plot([-1 3], [-1 3], 'k-', 'LineWidth', 2)
hold on
plot(y_oht, y_oht2, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'b')
xlabel('\Phi (m) 1 m Discretization')
ylabel('\Phi (m) 2 m Discretization')
ax.FontSize = 30;
