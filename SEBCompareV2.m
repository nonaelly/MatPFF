clear; close all; clc
addpath("Func\")

% Triangle element

%% Model 1
YourModel = 'Job-CZM-Tri.inp';  % Choose your model

sigma_c = 3.56; % MPa
G_c = 344; % J*m^-2
duP = 2.5e-3;
maxCMOD = 0.25; % mm

% lamda_cr = 0.001
[CMOD{1}, P{1}, u{1}, nodeEla{1}, elemEla{1}] = SEBFunctionTri(YourModel, 'bilinear', sigma_c, G_c, duP, maxCMOD, 0.001);

% lamda_cr = 0.04
[CMOD{2}, P{2}, u{2}] = SEBFunctionTri(YourModel, 'bilinear', sigma_c, G_c, duP, maxCMOD, 0.04);

% Exponential
[CMOD{3}, P{3}, u{3}] = SEBFunctionTri(YourModel, 'Exp', sigma_c, G_c, duP, maxCMOD, 0);

%% Plot
load Color_Config
figure
for i = 1 : 3
    if i < 4
        plot(CMOD{i}, - P{i},'MarkerEdgeColor', Color_Config{i}, 'LineWidth', 2)
    else
        plot(CMOD{i}(1 : 5 : end), - P{i}(1 : 5 : end),'MarkerEdgeColor', Color_Config{i}, ...
            'Marker', 'o', 'LineStyle', 'none', 'MarkerSize', 10, 'LineWidth', 1.5)
    end
    hold on
end

xlim([0 0.25])
lgd = legend('$${\rm{{\lambda}\ = 0.001 (Model 1)}}$$', '$${\rm{{\lambda}\ = 0.04 (Model 1)}}$$'...
    , '$${\rm{Exponential (Model 1)}}$$');
xtitle = '$${\mathrm{CMOD (mm)}}$$';
ytitle = '$${\mathrm{P (kN)}} $$';
setPlotV2(xtitle, ytitle, lgd)


%% Stress
dispV = reshape(u{1}(1 : size(nodeEla{1}, 1)*2, end), 2, [])';
% uy
figure
axis equal;
PlotContour(nodeEla{1},elemEla{1},dispV(:, 2),'uy',1);
axis off;

% ux
figure
axis equal;
PlotContour(nodeEla{1},elemEla{1},dispV(:, 1),'ux',1);
axis off;
