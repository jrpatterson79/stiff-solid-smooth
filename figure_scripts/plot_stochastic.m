% Oscillatory Flow Testing Fracture Characterization
% Numerical Modeling Analysis

% This code generates Figure 2 as seen in:
% Patterson, Jeremy R. & Cardiff, Michael (2023). Stiff, Smooth, and Solid?: Complex Fracture Hydraulic Hydraulics' Imprints on Oscillatory Hydraulic Testing. Submitted to Water Resources Research.

% Code developed by Jeremy Patterson
% Created June 2021; Updated Jan 2023

%% Clean Environment
close all; clear; clc

%% Load Files
file_dir = '/.../.../'; % Provide the directory location for the filename below

filename = 'stoch_exp_var_aniso.mat';
load([file_dir filename])

%% Stochastic Figures
idx = [find(syn_data(:,4)==r(4)) find(syn_data(:,4)==r(5)) find(syn_data(:,4)==r(2)) find(syn_data(:,4)==r(3))]; 

for k = 1:num_relz
    T(:,k) = [s_opt{k}(idx(:,4),1); s_opt{k}(idx(:,3),1)];
    S(:,k) = [s_opt{k}(idx(:,4),2); s_opt{k}(idx(:,3),2)];
    D(:,k) = [s_opt{k}(idx(:,4),3); s_opt{k}(idx(:,3),3)];
end

% Figure 7 - Stochastic results at 10 m
figure(7)
clf 
subplot(1,3,1)
ax = gca;
plot(P, exp(D(1:numel(P),:)), '-', 'Color', [0 0 0]+0.7, 'LineWidth', 0.5)
hold on
p1 = plot(P, exp(prctile(D(1:numel(P),:),2.5,2)), 'k', P, exp(prctile(D(1:numel(P),:),97.5,2)), 'k', 'LineWidth', 3);
p2 = plot(P, exp(mean(D(1:numel(P),:),2)), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 3);
p3 = plot([P(1) P(end)], [exp(lnD_mean) exp(lnD_mean)], '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 3);
ax.XScale = 'log';
ax.YScale = 'log';
xlim([P(1) P(end)])
xlabel('Period (s)')
ylabel('D (m^2/s)')
ax.FontSize = 30;
l = legend([p1(1), p2], '95% CI', 'Ensemble Mean', 'True Mean');
l.Location = 'NorthWest';
l.FontSize = 20;
text(700, 3e9, 'A', 'FontSize', 30, 'FontWeight', 'bold')

subplot(1,3,2)
ax = gca;
plot(P, exp(T(1:numel(P),:)), '-', 'Color', [0 0 0]+0.7, 'LineWidth', 0.5)
hold on
p1 = plot(P, exp(prctile(T(1:numel(P),:),2.5,2)), 'k', P, exp(prctile(T(1:numel(P),:),97.5,2)), 'k', 'LineWidth', 3);
p2 = plot(P, exp(mean(T(1:numel(P),:),2)), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 3);
p3 = plot([P(1) P(end)], [exp(lnT_mean) exp(lnT_mean)], '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 3);
axis([P(1) P(end) 1e-6 1e-4])
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('Period (s)')
ylabel('T (m^2/s)')
ax.FontSize = 30;
text(700, 8e-5, 'B', 'FontSize', 30, 'FontWeight', 'bold')

subplot(1,3,3)
ax = gca;
plot(P, exp(S(1:numel(P),:)), '-', 'Color', [0 0 0]+0.7, 'LineWidth', 0.5)
hold on
p1 = plot(P, exp(prctile(S(1:numel(P),:),2.5,2)), 'k', P, exp(prctile(S(1:numel(P),:),97.5,2)), 'k', 'LineWidth', 2);
p2 = plot(P, exp(mean(S(1:numel(P),:),2)), '-', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
p3 = plot([P(1) P(end)], [exp(lnS_mean) exp(lnS_mean)], '--', 'Color', [0 0.4470 0.7410], 'LineWidth', 3);
xlim([P(1) P(end)])
ax.XScale = 'log';
ax.YScale = 'log';
xlabel('Period (s)')
ylabel('S (-)')
ax.FontSize = 30;
text(700, 3e-6, 'C', 'FontSize', 30, 'FontWeight', 'bold')
set(gcf, 'Position', [0 0 2025 2025/3.5])