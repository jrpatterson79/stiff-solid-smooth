% Oscillatory Flow Testing Fracture Characterization
% Numerical Modeling Analysis

% This code generates Figure 2 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Stiff, Smooth, and Solid?: Complex Fracture Hydraulic Hydraulics' Imprints on Oscillatory Hydraulic Testing. Submitted to Water Resources Research.

% Code developed by Jeremy Patterson
% Created June 2021; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Load Heterogeneity Data
file_dir = '/.../.../'; % Provide the directory location to file1 and file2 immediately below

file1 = 'seed_50.mat';
file2 = 'seed_15.mat';

% Realization 1
load([file_dir file1])
aper_relz1 = log10(exp(ln_aper_true));
pen_dep_relz1 = pen_dep;
aper_mean_relz1 = log10(exp(aper_geo_mean));
y_relz1 = y; y_opt_relz1 = y_opt;
rmse_1 = sqrt(mean(y_relz1 - y_opt_relz1).^2);

% Realization 2
load([file_dir file2])
aper_relz2 = log10(exp(ln_aper_true));
pen_dep_relz2 = pen_dep;
aper_mean_relz2 = log10(exp(aper_geo_mean));
y_relz2 = y; y_opt_relz2 = y_opt;
rmse_2 = sqrt(mean(y_relz2 - y_opt_relz2).^2);

% Figure 2 - Aperture heterogeneity realizations
figure(2)
clf
tl = tiledlayout(2,2);

ax1 = nexttile(1);
p = pcolor(cgrid{1}, cgrid{2}, reshape(aper_relz1, numy, numx));
p.LineStyle = 'none';
hold on
plot(well_locs(1,1), well_locs(1,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', [0.7529 0 0])
% plot(well_locs(2:end,1), well_locs(2:end,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
axis([-50 50 -50 50])
axis square
ax1.XTick = [-50:25:50];
ax1.YTick = [-50:25:50];
xlabel('X distance (m)')
ylabel('Y distance (m)')
ax1.FontSize = 30;
ax1.YAxisLocation = 'right';
caxis([-5 -2]);
title('Realization 1', 'FontSize', 24, 'FontWeight', 'bold')
text(40, 45, 'A', 'FontSize', 30, 'FontWeight', 'bold')

ax3 = nexttile(3);
p = pcolor(cgrid{1}, cgrid{2}, reshape(aper_relz2, numy, numx));
p.LineStyle = 'none';
hold on
plot(well_locs(1,1), well_locs(1,2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', [0.7529 0 0])
axis([-50 50 -50 50])
axis square
ax3.XTick = [-50:25:50];
ax3.YTick = [-50:25:50];
xlabel('X distance (m)')
ylabel('Y distance (m)')
ax3.FontSize = 30;
ax3.YAxisLocation = 'right';
caxis([-5 -2])
title('Realization 2', 'FontSize', 24, 'FontWeight', 'bold')
text(40, 45, 'C', 'FontSize', 30, 'FontWeight', 'bold')

ax2 = nexttile(2);
plot(pen_dep_relz1, aper_mean_relz1, 'ko', 'MarkerFaceColor', [0.8500 0.3250 0.0980],...
    'MarkerSize', 10)
ylim([-3.55 -3.05])
axis square
ax2.XScale = 'log';
ax2.YAxisLocation = 'right';
xlabel('Distance from Pumping Well (m)')
ylabel('Mean log_{10}(Aperture [m])')
ax2.FontSize = 30;
text(650, -3.08, 'B', 'FontSize', 30, 'FontWeight', 'Bold')

ax4 = nexttile(4);
plot(pen_dep_relz2, aper_mean_relz2, 'ko', 'MarkerFaceColor', [0.8500 0.3250 0.0980],...
    'MarkerSize', 10)
ylim([-3.65 -3.45])
axis square
ax4.XScale = 'log';
ax4.YAxisLocation = 'right';
ax4.YTick = [-3.65:1e-1:-3.45];
xlabel('Distance from Pumping Well (m)')
ylabel('Mean log_{10}(Aperture [m])')
ax4.FontSize = 30;
text(650, -3.46, 'D', 'FontSize', 30, 'FontWeight', 'Bold')

c = colorbar;
caxis([-5 -2])
c.Label.String = 'log_{10}(Aperture [m])';
c.Layout.Tile = 'west';
tl.TileSpacing = 'loose';
tl.Padding = 'loose';
set(gcf, 'Position', [100 100 1100 1100/1.1])