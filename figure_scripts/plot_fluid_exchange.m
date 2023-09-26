% Oscillatory Flow Testing Fracture Characterization
% Numerical Modeling Analysis

% This code generates Figure 8 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Stiff, Smooth, and Solid?: Complex Fracture Hydraulic Hydraulics' Imprints on Oscillatory Hydraulic Testing. Submitted to Water Resources Research.

% Code developed by Jeremy Patterson
% Created June 2021; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Specify Directory
file_dir = '/.../.../'; % Provide the directory location to the files specified by file1, file2, and file3 below

%% Load Parameter Estimation Data
file1 = 'seed_50_3D.mat';
file2 = 'seed_15_3D.mat';
file3 = 'homog_3D.mat';

load([file_dir file1])
s_hat_relz1 = [s_opt(:,3) s_opt(:,1) s_opt(:,2)]; %[D,T,S]
y_relz1 = y_oht; y_opt_relz1 = [y_opt(:,1); y_opt(:,2)];
rmse_relz1 = sqrt(mean(y_opt_relz1 - y_relz1).^2);
ln_aper_relz1 = ln_aper_true;

load([file_dir file2])
s_hat_relz2 = [s_opt(:,3) s_opt(:,1) s_opt(:,2)]; %[D,T,S]
y_relz2 = y_oht; y_opt_relz2 = [y_opt(:,1); y_opt(:,2)];
rmse_relz2 = sqrt(mean(y_opt_relz2 - y_relz2).^2);
lnK_relz2 = lnK_true;

load([file_dir file3])
s_hat_homog = [s_opt(:,3) s_opt(:,1) s_opt(:,2)]; %[D,T,S]
y_homog = y_oht; y_opt_homog = [y_opt(:,1); y_opt(:,2)];
rmse_homog = sqrt(mean(y_opt_homog - y_homog).^2);

% Calculate mean parameter values
lnT_mean = log(((rho*g)/(12 * mu)) .* exp(ln_aper_mean).^3);
lnS_mean = log(Ss_eff * exp(ln_aper_mean));
lnD_mean = lnT_mean - lnS_mean;

% Calculate correlation coefficient for each parameter trend with period
% and distance as well as p-value testing the hypothesis that the
% correlation is different from correlation of 0. 
% First column is parameter correlation with period and p-values
% Second column is parameter correlation with distance and p-values
% Row 1 is diffusivity, row 2 is transmissivity, row 3 is storativity
rho = zeros(3,2); pval = zeros(3,2);
for i = 1:3
    % Realization 1
    [rho_relz1(i,1), pval_relz1(i,1)] = corr(log(syn_data(:,1)), s_hat_relz1(:,i));
    [rho_relz1(i,2), pval_relz1(i,2)] = corr(log(syn_data(:,4)), s_hat_relz1(:,i));
    % Realization 2
    [rho_relz2(i,1), pval_relz2(i,1)] = corr(log(syn_data(:,1)), s_hat_relz2(:,i));
    [rho_relz2(i,2), pval_relz2(i,2)] = corr(log(syn_data(:,4)), s_hat_relz2(:,i));
end

%% Figures
% Uncomment the lines below to plot volume slices that show the modeling domain
% lbls = {'A1', 'B3', 'B4', 'B1', 'B2'};
% figure
% clf
% ax = gca;
% sx = slice(cgrid{1}, cgrid{2}, cgrid{3}, reshape(log10(exp(ln_aper_relz1)), numy, numx, numz), 15, [], []);
% hold on
% sy = slice(cgrid{1}, cgrid{2}, cgrid{3}, reshape(log10(exp(ln_aper_relz1)), numy, numx, numz), [], 15, []);
% sz = slice(cgrid{1}, cgrid{2}, cgrid{3}, reshape(log10(exp(ln_aper_relz1)), numy, numx, numz), [], [], node_list(1));
% sx.FaceAlpha = 0.9;
% sy.FaceAlpha = 0.9;
% sx.FaceColor = 'flat';
% sy.FaceColor = 'flat';
% sz.FaceColor = 'flat';
% plot3(well_locs(:,1), well_locs(:,2), well_locs(:,3), 'ko',...
%       'MarkerFaceColor', [192/255 0 0], 'MarkerSize', 14)
% text(well_locs(:,1), well_locs(:,2), well_locs(:,3), lbls, 'FontSize', 24,...
%      'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')
% axis([-40 40 -40 40 node_list(1) zmax])
% ax.XTick = [-25:10:25];
% ax.YTick = [-25:10:25];
% ax.ZTick = [0:0.1:0.2];
% view(-35, 60)
% c1 = colorbar;
% caxis([-2 -5])
% % c1.Ticks = [-15:3:-3];
% c1.Label.String = 'log_{10}(K [m/s])';
% c1.FontSize = 24;
% xlabel('X (m)')
% ylabel('Y (m)')
% zlabel('Z (m)')
% ax.FontSize = 28;
% set(gcf, 'Position', [100 100 975 975/1.5])


% Figure 6 - Panels G, H, I (See plot_figure_6.m to reproduce full figure)  
figure
clf
subplot(2,3,[1:2,4:5])
ax = gca;
scatter(syn_data(:,1), exp(s_hat_homog(:,1)), 150, syn_data(:,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(lnD_mean) exp(lnD_mean)], 'k--', 'LineWidth', 3)
axis([1e1 1e3 1e-4 1e1])
ax.XScale = 'log';
ax.YScale = 'log';
grid on
xlabel('Oscillation Period (s)')
ylabel('D_{est} (m^2/s)')
ax.FontSize = 30;
c = colorbar;
caxis([4 16])
c.Label.String = 'Radial Distance (m)';
c.FontSize = 24;
text(10.75, 7e0, 'A', 'FontSize', 30, 'FontWeight', 'bold')

subplot(2,3,3)
ax = gca;
scatter(syn_data(:,1), exp(s_hat_homog(:,2)), 150, syn_data(:,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(lnT_mean) exp(lnT_mean)], 'k--', 'LineWidth', 3)
ax.XScale = 'log';
ax.YScale = 'log';
ax.YTick = [1e-9; 1e-8; 1e-7; 1e-6; 1e-5];
axis([1e1 1e3 1e-9 1e-5])
grid on
caxis([4 16])
ylabel('T_{est} (m^2/s)')
ax.YAxisLocation = 'right';
ax.FontSize = 30;
text(10.75, 6e-6, 'B', 'FontSize', 30, 'FontWeight', 'bold')

subplot(2,3,6)
ax = gca;
scatter(syn_data(:,1), exp(s_hat_homog(:,3)), 150, syn_data(:,4),...
        'filled', 'MarkerEdgeColor', 'k')
hold on
plot([P(1) P(end)], [exp(lnS_mean) exp(lnS_mean)], 'k--', 'LineWidth', 3)
axis([1e1 1e3 1e-7 1e-4])
ax.XScale = 'log';
ax.YScale = 'log';
grid on
caxis([4 16])
xlabel('Oscillation Period (s)')    
ylabel('S_{est} (-)')
ax.YAxisLocation = 'right';
ax.FontSize = 30;
text(10.75, 6e-5, 'C', 'FontSize', 30, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 2025 2025/1.8])

% Figure 8 - Heterogeneity and fluid exchange
% Selects a subset of inter-well spacings for figure clarity. Replace (idx,#) with (:,#) in the figure below to show all distances
idx = find(syn_data(:,4) == 4 | syn_data(:,4) == 8 | syn_data(:,4) == 16);
figure(8)
clf
subplot(2,3,[1:2,4:5])
% Hydraulic Diffusivity
ax1 = gca;
scatter(syn_data(idx,1), exp(s_hat_relz1(idx,1)), 200, syn_data(idx,4),...
        's', 'filled', 'MarkerEdgeColor', 'k')
hold on
scatter(syn_data(idx,1), exp(s_hat_relz2(idx,1)), 200, syn_data(idx,4),...
        '>', 'filled', 'MarkerEdgeColor', 'k')
plot([P(1) P(end)], [exp(lnD_mean) exp(lnD_mean)], 'k--', 'LineWidth', 3)
ax1.XScale = 'log';
ax1.YScale = 'log';
axis([1e1 1e3 1e-4 1e1])
grid on
grid minor
ax1.MinorGridAlpha = 0.7;
caxis([4 10])
xlabel('Oscillation Period (s)')
ylabel('D_{est} (m^2/s)')
ax1.FontSize = 30;
[l, obj] = legend('Realization 1', 'Realization 2');
objl = findobj(obj, 'type', 'patch');
set(objl, 'Markersize', 12)
l.Box = 'off';
l.FontSize = 20;
l.Location = 'northwest';
c = colorbar;
caxis([4 16])
c.Label.String = 'Inter-well Spacing (m)';
c.FontSize = 24;
text(800, 7e0, 'A', 'FontSize', 30, 'FontWeight', 'bold')

% Transmissivity
subplot(2,3,3)
ax2 = gca;
scatter(syn_data(idx,1), exp(s_hat_relz1(idx,2)), 200, syn_data(idx,4),...
        's', 'filled', 'MarkerEdgeColor', 'k')
hold on
scatter(syn_data(idx,1), exp(s_hat_relz2(idx,2)), 200, syn_data(idx,4),...
        '>', 'filled', 'MarkerEdgeColor', 'k')
plot([P(1) P(end)], [exp(lnT_mean) exp(lnT_mean)], 'k--', 'LineWidth', 3)
axis([1e1 1e3 1e-9 1e-4])
ax2.XScale = 'log';
ax2.YScale = 'log';
ax.YTick = [1e-9; 1e-8; 1e-7; 1e-6; 1e-5; 1e-4];
grid on
ax2.MinorGridAlpha = 0.7;
caxis([4 16])
ylabel('T_{est} (m^2/s)')
ax2.YAxisLocation = 'right';
ax2.FontSize = 30;
text(725, 6e-5, 'B', 'FontSize', 30, 'FontWeight', 'bold')

% Storativity
subplot(2,3,6)
ax3 = gca;
scatter(syn_data(idx,1), exp(s_hat_relz1(idx,3)), 200, syn_data(idx,4),...
        's', 'filled', 'MarkerEdgeColor', 'k')
hold on
scatter(syn_data(idx,1), exp(s_hat_relz2(idx,3)), 200, syn_data(idx,4),...
        '>', 'filled', 'MarkerEdgeColor', 'k')
plot([P(1) P(end)], [exp(lnS_mean) exp(lnS_mean)], 'k--', 'LineWidth', 3)
axis([1e1 1e3 1e-8 1e-4])
ax3.XScale = 'log';
ax3.YScale = 'log';
grid on
ax3.MinorGridAlpha = 0.7;
caxis([4 16])
xlabel('Oscillation Period (s)')    
ylabel('S_{est} (-)')
ax3.YAxisLocation = 'right';
ax3.FontSize = 30;
text(725, 7e-5, 'C', 'FontSize', 30, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 2025 2025/1.8])