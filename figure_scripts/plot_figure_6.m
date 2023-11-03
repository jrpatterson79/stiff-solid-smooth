% Oscillatory Flow Testing Fracture Characterization
% Numerical Modeling Analysis

% This code generates Figure 6 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Stiff, Smooth, and Solid?: Complex Fracture Hydraulic Hydraulics' Imprints on Oscillatory Hydraulic Testing. Submitted to Water Resources Research.

% Code developed by Jeremy Patterson
% Created Nov 2022; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Load Field Data

field_dir = '/.../.../'; % Provide the directory location to '2019_pest_results.mat'
load([field_dir '2019_pest_results.mat'])

% Find field tests with pumping period >= 10 s
field_idx = find(s_unc_19(:,1) < 1 & test_list(:,1) > 4);
field_list = test_list;
s_hat_field = [s_opt_19(:,3) s_opt_19(:,1) s_opt_19(:,2)];
clear test_list

%% Load Heterogeneous Modeling Results
% Provide the directory to .mat files with heterogeneity results 'seed_50.mat' and 'seed_15.mat'
heterog_dir = '/.../.../';

% Realization 1
load([heterog_dir 'seed_50.mat'])
heterog_list = syn_data;
lnT_true_relz1 = lnT_true;
s_hat_relz1 = [s_hat(:,1)-s_hat(:,2) s_hat(:,1) s_hat(:,2)];

% Realization 2
load([heterog_dir 'seed_15.mat'])
s_hat_relz2 = [s_hat(:,1)-s_hat(:,2) s_hat(:,1) s_hat(:,2)];

%% Load Fracture-Host Rock Interaction Modeling Results
lky_dir = '/.../.../'; % Provide the directory location to the fluid exchange results 'homog_3D.mat'
load([lky_dir 'homog_3D.mat'])
s_hat_lky = [s_opt(:,3) s_opt(:,1) s_opt(:,2)];
lky_list = syn_data;

%% Load Hydromechanical Modeling Results
% mech_dir = '/.../.../'; % Provide the directory location the the hydromechanical results 'poromech_pest_results.mat'
load([mech_dir 'hydromechanical_pest_results.mat'])

s_hat_mech = [s_hat_imperm(:,3) s_hat_imperm(:,1) s_hat_imperm(:,2)];
mech_list = test_list;

%% Table 5
% Calculate correlation coefficient for each parameter trend with period
% and distance as well as p-value testing (two-sided) the hypothesis that the
% correlation is different from correlation of 0. 
% First column is parameter correlation with period and p-values
% Second column is parameter correlation with distance and p-values
% Row 1 is diffusivity, row 2 is transmissivity, row 3 is storativity

rho_field = zeros(3,2); rho_relz1 = zeros(3,2); rho_relz2 = zeros(3,2); rho_lky = zeros(3,2); rho_mech = zeros(3,2);
pval_field = zeros(3,2); pval_relz1 = zeros(3,2); pval_relz2 = zeros(3,2); pval_lky = zeros(3,2); pval_mech = zeros(3,2); 
for i = 1:3
    [rho_field(i,1), ~] = corr(log(field_list(field_idx,1)), s_hat_field(field_idx,i));
    [rho_field(i,2), ~] = corr(log(field_list(field_idx,4)), s_hat_field(field_idx,i));
    [rho_relz1(i,1), ~] = corr(log(heterog_list(:,1)), s_hat_relz1(:,i));
    [rho_relz1(i,2), ~] = corr(log(heterog_list(:,4)), s_hat_relz1(:,i));
    [rho_relz2(i,1), ~] = corr(log(heterog_list(:,1)), s_hat_relz2(:,i));
    [rho_relz2(i,2), ~] = corr(log(heterog_list(:,4)), s_hat_relz2(:,i));
    [rho_lky(i,1), ~] = corr(log(lky_list(:,1)), s_hat_lky(:,i));
    [rho_lky(i,2), ~] = corr(log(lky_list(:,4)), s_hat_lky(:,i));
    [rho_mech(i,1), ~] = corr(log(mech_list(:,1)), s_hat_mech(:,i));
    [rho_mech(i,2), ~] = corr(log(mech_list(:,4)), s_hat_mech(:,i));
end

%% Figures
% Figure 6 - Multi-panel plot of parameter estimation results from all
% models / explored mechanisms
figure(6)
clf
tl = tiledlayout(3,4);
% Field Diffusivity
ax1 = nexttile(1);
scatter(field_list(field_idx,1), exp(s_hat_field(field_idx,1)), 125, field_list(field_idx,4),...
        'filled', 'MarkerEdgeColor', 'k')
grid on
ax1.MinorGridAlpha = 0.7;
axis([min(field_list(field_idx,1)) max(field_list(field_idx,1)) 1e-2 1e2])
caxis([4 16])
ax1.YTick = [1e-2 1e0 1e2];
ax1.XScale = 'log';
ax1.YScale = 'log';
ylabel('Diffusivity (m^2/s)')
ax1.FontSize = 30;
text(225, 65, 'A', 'FontSize', 30, 'FontWeight', 'bold')
% title('Field Results', 'FontSize', 24)

% Field Transmissivity
ax5 = nexttile(5);
scatter(field_list(field_idx,1), exp(s_hat_field(field_idx,2)), 125, field_list(field_idx,4),...
        'filled', 'MarkerEdgeColor', 'k')
grid on
ax5.MinorGridAlpha = 0.7;
axis([min(field_list(field_idx,1)) max(field_list(field_idx,1)) 1e-7 1e-3])
ax5.YTick = [1e-7;1e-5;1e-3];
ax5.XScale = 'log';
ax5.YScale = 'log';
caxis([4 16])
ylabel('Transmissivity (m^2/s)')
ax5.FontSize = 30;
text(225, 6.5e-4, 'B', 'FontSize', 30, 'FontWeight', 'bold')

% Field Storativity
ax9 = nexttile(9);
scatter(field_list(field_idx,1), exp(s_hat_field(field_idx,3)), 125, field_list(field_idx,4),...
        'filled', 'MarkerEdgeColor', 'k')
grid on
ax9.MinorGridAlpha = 0.7;
axis([min(field_list(field_idx,1)) max(field_list(field_idx,1)) 1e-7 1e-3])
ax9.YTick = [1e-7; 1e-5; 1e-3];
ax9.XScale = 'log';
ax9.YScale = 'log';
caxis([4 16])
ylabel('Storativity (-)')
ax9.FontSize = 30;
text(225, 6.5e-4, 'C', 'FontSize', 30, 'FontWeight', 'bold')

% Heterogeneity Diffusivity
ax2 = nexttile(2);
scatter(heterog_list(:,1), exp(s_hat_relz1(:,1)), 125, heterog_list(:,4),...
        's', 'filled', 'MarkerEdgeColor', 'k')
hold on
scatter(heterog_list(:,1), exp(s_hat_relz2(:,1)), 125, heterog_list(:,4),...
        '>', 'filled', 'MarkerEdgeColor', 'k')
grid on
ax2.MinorGridAlpha = 0.7;
ylim([1e0 1e4])
caxis([4 16])
ax2.YTick = [1e0 1e2 1e4];
ax2.XScale = 'log';
ax2.YScale = 'log';
ax2.FontSize = 30;
[l, obj] = legend('Realization 1', 'Realization 2');
objl = findobj(obj, 'type', 'patch');
set(objl, 'Markersize', 12)
l.Box = 'off';
l.FontSize = 20;
l.Location = 'northwest';
text(650, 6.5e3, 'D', 'FontSize', 30, 'FontWeight', 'bold')
% title('Aperture Heterogeneity', 'FontSize', 24)

% Heterogeneity Transmissivity
ax6 = nexttile(6);
scatter(heterog_list(:,1), exp(s_hat_relz1(:,2)), 125, heterog_list(:,4),...
        's', 'filled', 'MarkerEdgeColor', 'k')
hold on
scatter(heterog_list(:,1), exp(s_hat_relz2(:,2)), 125, heterog_list(:,4),...
        '>', 'filled', 'MarkerEdgeColor', 'k')
plot([P(1) P(end)], [exp(lnT_mean) exp(lnT_mean)], 'k--', 'LineWidth', 3)
grid on
ax6.MinorGridAlpha = 0.7;
ylim([1e-7 1e-3])
ax6.YTick = [1e-7; 1e-5; 1e-3];
caxis([4 16])
ax6.XScale = 'log';
ax6.YScale = 'log';
ax6.FontSize = 30;
text(650, 6.5e-4, 'E', 'FontSize', 30,'FontWeight', 'bold')

% Heterogeneity Storativity
ax10 = nexttile(10);
scatter(heterog_list(:,1), exp(s_hat_relz1(:,3)), 125, heterog_list(:,4),...
        's', 'filled', 'MarkerEdgeColor', 'k')
hold on
scatter(heterog_list(:,1), exp(s_hat_relz2(:,3)), 125, heterog_list(:,4),...
        '>', 'filled', 'MarkerEdgeColor', 'k')
plot([P(1) P(end)], [exp(lnS_mean) exp(lnS_mean)], 'k--', 'LineWidth', 3)
grid on
ax10.MinorGridAlpha = 0.7;
axis([1e1 1e3 1e-9 1e-5])
ax10.YTick = [1e-9;1e-7;1e-5];
caxis([4 16])
ax10.XScale = 'log';
ax10.YScale = 'log';
ax10.FontSize = 30;
text(650, 6.5e-6, 'F', 'FontSize', 30,'FontWeight', 'bold')

% Leaky Diffusivity
ax3 = nexttile(3);
scatter(lky_list(:,1), exp(s_hat_lky(:,1)), 125, lky_list(:,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(lnD_mean) exp(lnD_mean)], 'k--', 'LineWidth', 3)
grid on
ax3.MinorGridAlpha = 0.7;
axis([1e1 1e3 1e-3 1e1])
ax3.YTick = [1e-3; 1e-1; 1e1];
caxis([4 16])
ax3.XScale = 'log';
ax3.YScale = 'log';
ax3.FontSize = 30;
text(650, 6.5, 'G', 'FontSize', 30, 'FontWeight', 'bold')
% title('Fracture-Host Rock Fluid Exchange', 'FontSize', 24)

% Leaky Transmissivity
ax7 = nexttile(7);
scatter(lky_list(:,1), exp(s_hat_lky(:,2)), 125, lky_list(:,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(lnT_mean) exp(lnT_mean)], 'k--', 'LineWidth', 3)
grid on
ax7.MinorGridAlpha = 0.7;
axis([1e1 1e3 1e-8 1e-4])
caxis([4 16])
ax7.YTick = [1e-8; 1e-6; 1e-4];
ax7.XScale = 'log';
ax7.YScale = 'log';
ax7.FontSize = 30;
text(650, 6.5e-5, 'H', 'FontSize', 30,'FontWeight', 'bold')

% Leaky Storativity
ax11 = nexttile(11);
scatter(lky_list(:,1), exp(s_hat_lky(:,3)), 125, lky_list(:,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(lnS_mean) exp(lnS_mean)], 'k--', 'LineWidth', 3)
grid on
ax11.MinorGridAlpha = 0.7;
axis([1e1 1e3 1e-8 1e-4])
caxis([4 16])
ax11.YTick = [1e-8; 1e-6; 1e-4];
ax11.XScale = 'log';
ax11.YScale = 'log';
ax11.FontSize = 30;
text(650, 6.5e-5, 'I', 'FontSize', 30,'FontWeight', 'bold')

% Mechanical Diffusivity
ax4 = nexttile(4);
scatter(mech_list(:,1), exp(s_hat_mech(:,1)), 125, mech_list(:,4),...
        'filled', 'MarkerEdgeColor', 'k')
grid on
ax4.MinorGridAlpha = 0.7;
axis([1e1 1e3 1e-3 1e1])
ax4.YTick = [1e-3; 1e-1; 1e1];
caxis([4 16])
ax4.XScale ='log';
ax4.YScale = 'log';
ax4.FontSize = 30;
text(650, 6.5, 'J', 'FontSize', 30,'FontWeight', 'bold')
% title('Fracture Hydromechanics', 'FontSize', 24)

% Mechanical Transmissivity
ax8 = nexttile;
scatter(mech_list(:,1), exp(s_hat_mech(:,2)), 125, mech_list(:,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on 
plot([P(1) P(end)], [exp(lnT_mean) exp(lnT_mean)], 'k--', 'LineWidth', 3)
grid on
ax8.MinorGridAlpha = 0.7;
axis([1e1 1e3 1e-8 1e-4])
caxis([4 16])
ax8.YTick = [1e-8; 1e-6; 1e-4];
ax8.XScale ='log';
ax8.YScale = 'log';
ax8.FontSize = 30;
text(650, 6.5e-5, 'K', 'FontSize', 30,'FontWeight', 'bold')

% Mechanical Storativity
ax12 = nexttile(12);
scatter(mech_list(:,1), exp(s_hat_mech(:,3)), 125, mech_list(:,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(lnS_mean) exp(lnS_mean)], 'k--', 'LineWidth', 3)
grid on
ax12.MinorGridAlpha = 0.7;
axis([1e1 1e3 1e-8 1e-4])
caxis([4 16])
ax12.YTick = [1e-8; 1e-6; 1e-4];
ax12.XScale ='log';
ax12.YScale = 'log';
ax12.FontSize = 30;
text(650, 6.5e-5, 'L', 'FontSize', 30,'FontWeight', 'bold')

xlabel(tl, 'Oscillation Period (s)', 'FontSize', 30)
tl.TileSpacing = 'loose';
tl.Padding = 'loose';
c = colorbar;
c.Layout.Tile = 'east';
c.Label.String = 'Inter-well Spacing (m)';
c.FontSize = 26;
set(gcf, 'Position', [0 0 2500 2500/1.1])
