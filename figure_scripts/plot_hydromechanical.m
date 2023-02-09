% Oscillatory Flow Testing Fracture Characterization
% Numerical Modeling Analysis

% This code generates Figure 5, Figure 9, Figure 10, and Figure 11 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Stiff, Solid, and Smooth?: Complex Fracture Hydraulic Hydraulics Revealed through Oscillatory Flow Interference Testing. Submitted to Water Resources Research

% Code developed by Jeremy Patterson
% Created August 2021; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Specify Directory
addpath('/.../.../') % Provide the directory location to the subdirectory func_lib which contains all of the necessary function files 

%% Import and Parse Data
mech_dir = '/.../.../'; % Specify the directory location to the hydromechanical mat file indicated on the line below
load([mech_dir 'hydromechanical_pest_results.mat'])

%% Figures
% Figure 5 - Pressure and displacement time series
p_idx = data.P_s_ == P(6);
figure(5)
clf
yyaxis left
ax = gca;
plot(data.time_s_(p_idx), h_obs(p_idx,4), '-', 'LineWidth', 3)
grid on
xlabel('Time (s)')
ylabel('Head Change (m)')
axis([0 P(6)*10 -0.15 0.15])
ax.FontSize = 30;

yyaxis right
plot(data.time_s_(p_idx), strain_obs(p_idx,4)*1e6, '-.', 'LineWidth', 3)
ylim([-0.52 0.52])
ax.YTick = [-0.5:0.25:0.5];
axis([0 P(6)*10 -0.55 0.55])
xlabel('Time (s)')
ylabel('Fracture Displacement (\mum)')
ax.FontSize = 32;
set(gcf, 'Position', [0 0 2025 2025/3.1])

% Figure 9 - Pressure / displacement amplitude subplot (Perm vs Imperm)
figure(9)
tl = tiledlayout(2,2);
% Impermeable Pressure Amplitude
ax1 = nexttile(1);
scatter(test_list(:,1), abs(press_phasor_imperm), 125, test_list(:,4), 'filled',...
        'MarkerEdgeColor', 'k')
ax1.XScale = 'log';
ax1.YScale = 'log';
grid on
axis([0 P(end) 2e-2 1e0])
caxis([4 16])
ylabel('Head Amplitude (m)')
ax1.FontSize = 30;
text(10.75, 0.75, 'A', 'FontSize', 30, 'FontWeight', 'bold')
title('Impermeable Host Rock', 'FontSize', 24)
axis square

% Permeable Pressure Amplitude
ax2 = nexttile(2);
scatter(test_list(:,1), abs(press_phasor_perm), 125, test_list(:,4), 'filled',...
        'MarkerEdgeColor', 'k')
ax2.XScale = 'log';
ax2.YScale = 'log';
grid on
axis([0 P(end) 2e-2 1e0])
caxis([4 16])
ylabel('Head Amplitude (m)')
ax2.FontSize = 30;
title('Permeable Host Rock', 'FontSize', 24)
text(10.75, 0.75, 'C', 'FontSize', 30, 'FontWeight', 'bold')
axis square

% Impermeable Strain Amplitude
ax3 = nexttile(3);
scatter(test_list(:,1), abs(strain_phasor_imperm)*1e6, 125, test_list(:,4), 'filled',...
        'MarkerEdgeColor', 'k')
ax3.XScale = 'log';
ax3.YScale = 'log';
grid on
axis([0 P(end) 4e-2 1e1])
caxis([4 16])
ylabel('Displacement Amplitude (\mum)')
ax3.FontSize = 30;
text(10.75, 7.5, 'B', 'FontSize', 30, 'FontWeight', 'bold')
axis square

% Permeable Strain Amplitude
ax4 = nexttile(4);
scatter(test_list(:,1), abs(strain_phasor_perm)*1e6, 125, test_list(:,4), 'filled',...
        'MarkerEdgeColor', 'k')
ax4.XScale = 'log';
ax4.YScale = 'log';
grid on
axis([0 P(end) 4e-2 1e1])
caxis([4 16])
ylabel('Displacement Amplitude (\mum)')
ax4.FontSize = 30;
text(10.75, 7.5, 'D', 'FontSize', 30, 'FontWeight', 'bold')
axis square

xlabel(tl, 'Oscillation Period (s)', 'FontSize', 30)
tl.TileSpacing = 'loose';
tl.Padding = 'loose';
c = colorbar;
c.Layout.Tile = 'east';
caxis([4 16])
c.Label.String = 'Inter-well Spacing (m)';
c.FontSize = 24;
set(gcf, 'Position', [0 0 1100 1100/1.1])

% Figure 10 - Impermeable vs Permeable parameter estimates
figure(10)
clf
subplot(2,3,[1:2,4:5])
ax = gca;
scatter(test_list(:,1), exp(s_hat_perm(:,3)), 125, test_list(:,4),...
        'd', 'filled', 'MarkerEdgeColor', 'k')
grid on
axis([1e1 1e3 1e-3 2e0])
ax.XScale ='log';
ax.YScale = 'log';
xlabel('Oscillation Period (s)')
ylabel('D_{est} (m^2/s)')
ax.FontSize = 30;
c = colorbar;
caxis([4 16])
c.Label.String = 'Inter-well Spacing (m)';
c.FontSize = 24;
text(800, 1.7, 'A', 'FontSize', 30, 'FontWeight', 'bold')

subplot(2,3,3)
ax = gca;
scatter(test_list(:,1), exp(s_hat_perm(:,1)), 125, test_list(:,4),...
        'd', 'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(lnT) exp(lnT)], 'k--', 'LineWidth', 3)
grid on
axis([1e1 1e3 1e-8 1e-4])
ax.XScale = 'log';
ax.YScale = 'log';
caxis([4 16])
ylabel('T_{est} (m^2/s)')
ax.YAxisLocation = 'right';
ax.FontSize = 30;
text(650, 6.5e-5, 'B', 'FontSize', 30, 'FontWeight', 'bold')

subplot(2,3,6)
ax = gca;
scatter(test_list(:,1), exp(s_hat_perm(:,2)), 125, test_list(:,4),...
        'd', 'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(-17.4831) exp(-17.4831)], 'k--', 'LineWidth', 3)
grid on
ax.XScale = 'log';
ax.YScale = 'log';
axis([1e1 1e3 1e-8 1e-4])
caxis([4 16])
xlabel('Oscillation Period (s)')    
ylabel('S_{est} (-)')
ax.YAxisLocation = 'right';
ax.FontSize = 30;
text(650, 6.5e-5, 'C', 'FontSize', 30, 'FontWeight', 'bold')
set(gcf, 'Position', [100 100 2025 2025/1.8])

% Figure 11 - Effective T for multiple k
r_idx = test_list(:,4) == 4;
lbls = {'K = 1e-15 m/s', 'K = 3e-8 m/s', 'K = 3e-6 m/s'};
col = {[0.9290 0.6940 0.1250], [0.4660 0.6740 0.1880], [0.6350 0.0780 0.1840]};

figure(11)
clf
ax = gca;
plot(test_list(r_idx,1), exp(s_hat_imperm(r_idx,1)), 'ko', 'MarkerFaceColor', col{1}, 'MarkerSize', 12)
hold on
plot(test_list(r_idx,1), exp(s_hat_perm(r_idx,1)), 'ks', 'MarkerFaceColor', col{2}, 'MarkerSize', 12)
plot(test_list(r_idx,1), exp(s_hat_13(r_idx,1)), 'k>', 'MarkerFaceColor', col{3}, 'MarkerSize', 12)
grid on
ylim([1e-8 1e-6])
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('Oscillation Period (s)')
ylabel('ln(T_{est} [m^2/s])')
legend(lbls)
ax.FontSize = 30;
set(gcf, 'Position', [100 100 975 975/1.3333])