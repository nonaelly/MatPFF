clear; close all; clc
addpath("Func\")

%% Model 1
YourModel = 'Job-CZM-mesh.inp';  % Choose your model
% YourModel = 'Job-CZM-mesh2.inp';  % Choose your model

sigma_c = 3.56; % MPa
G_c = 344; % J*m^-2
duP = 2.5e-3;
maxCMOD = 0.25; % mm

% lamda_cr = 0.001
[CMOD{1}, P{1}, u{1}, nodeEla{1}, elemEla{1}] = SEBFunction(YourModel, 'bilinear', sigma_c, G_c, duP, maxCMOD, 0.001);
% dispV = reshape(u{1}(1 : size(nodeEla{1}, 1)*2, end), 2, [])';
% vtkName = ['lamda=1','.vtk'];
% vtkwrite(vtkName, 'UNSTRUCTURED_GRID', nodeEla{1}(:,1), nodeEla{1}(:,2), nodeEla{1}(:,1)*0, ...
%     'cells', elemEla{1}, 'cell_types', 9, 'vectors', 'u', dispV(:, 1)', dispV(:, 2)', 0*dispV(:, 1)');

% lamda_cr = 0.04
[CMOD{2}, P{2}, u{2}] = SEBFunction(YourModel, 'bilinear', sigma_c, G_c, duP, maxCMOD, 0.04);

% Exponential
[CMOD{3}, P{3}, u{3}] = SEBFunction(YourModel, 'Exp', sigma_c, G_c, duP, maxCMOD, 0);

%% Model 2
% YourModel = 'Job-CZM-mesh.inp';  % Choose your model
YourModel = 'Job-CZM-mesh2.inp';  % Choose your model

sigma_c = 3.56; % MPa
G_c = 344; % J*m^-2
duP = 2.5e-3;
maxCMOD = 0.25; % mm

% lamda_cr = 0.001
[CMOD{4}, P{4}, u{4}, nodeEla{2}, elemEla{2}] = SEBFunction(YourModel, 'bilinear', sigma_c, G_c, duP, maxCMOD, 0.001);

% lamda_cr = 0.04
[CMOD{5}, P{5}, u{5}] = SEBFunction(YourModel, 'bilinear', sigma_c, G_c, duP, maxCMOD, 0.04);

% Exponential
[CMOD{6}, P{6}, u{6}] = SEBFunction(YourModel, 'Exp', sigma_c, G_c, duP, maxCMOD, 0);

%% Plot
load Color_Config
figure
% for i = 1 : 6
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
ylim([0 9])
lgd = legend('$${\rm{{\lambda}\ = 0.001 (Model 1)}}$$', '$${\rm{{\lambda}\ = 0.04 (Model 1)}}$$'...
    , '$${\rm{Exponential (Model 1)}}$$', '$${\rm{{\lambda}\ = 0.001 (Model 2)}}$$' ...
    , '$${\rm{{\lambda}\ = 0.04 (Model 2)}}$$', '$${\rm{Exponential (Model 2)}}$$');
% lgd = legend('$${\rm{{\lambda}\ = 0.001 (Model 1)}}$$', '$${\rm{{\lambda}\ = 0.04 (Model 1)}}$$'...
%     , '$${\rm{Exponential (Model 1)}}$$');
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


% dispV = reshape(u{4}(1 : size(nodeEla{2}, 1)*2, end), 2, [])';
% % uy
% figure
% axis equal;
% PlotContour(nodeEla{2},elemEla{2},dispV(:, 2),'uy',1);
% axis off;
% 
% % ux
% figure
% axis equal;
% PlotContour(nodeEla{2},elemEla{2},dispV(:, 1),'ux',1);
% axis off;
