clear; close all; clc
addpath("Func\")

%% Model 1
% YourModel = 'Job-bulk2.inp';  % Choose your model
YourModel = 'Job-CZM-Bulk.inp';  % Choose your model

sigma_c = 3.56; % MPa
G_c = 344; % J*m^-2
duP = 2.5e-3;
maxCMOD = 0.25; % mm

% lamda_cr = 0.001
%%  ***  Read Abaqus Mesh  ***
parts = loadinp(YourModel);

% Copy the half model
nodes = [(1 : size(parts.Nodes, 1))', parts.Nodes];

elems = [];
for i = 1 : size(parts.Elements, 1)
    [elems] = [elems; i, 1, parts.Elements(i).Connectivity];
end

insertionRegions = [0, 0, 0, 19 - 1e-5];
[nodes, newElements, cohesiveElements] = insertCohesiveElements(nodes, elems, insertionRegions);
% Pre-define
elems = newElements;

insertionRegions = [-1e-5, 1e-5, 19 - 1e-5, 100 + 1e-5];
[nodes, newElements, cohesiveElements] = insertCohesiveElements(nodes, elems, insertionRegions);
elems = [newElements; cohesiveElements];

nodeX = nodes(:, 2);
nodeY = nodes(:, 3);

numElements = size(elems(elems(:, 2) == 1), 1);

X = zeros(numElements, 4); % 每个单元4个顶点
Y = zeros(numElements, 4);

for i = 1 : numElements
    elementNodes = elems(i, 3:6);  % 假设每个单元是四边形
    
    X(i, :) = nodeX(elementNodes);
    Y(i, :) = nodeY(elementNodes);
end

figure;
patch(X', Y', 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'k');
hold on

X = zeros(size(elems, 1) - numElements, 4); % 每个单元4个顶点
Y = zeros(size(elems, 1) - numElements, 4);

for i = numElements + 1 : size(elems, 1)
    elementNodes = elems(i, 3:6);  % 假设每个单元是四边形
    
    X(i, :) = nodeX(elementNodes);
    Y(i, :) = nodeY(elementNodes);
end

patch(X', Y', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'g');

axis equal;
title('Finite Element Mesh');
xlabel('X');
ylabel('Y');


title('single-edge notched beam (SE(B)) test');
xlabel('X');
ylabel('Y');
axis equal;

sigma_c = 3.56; % MPa
G_c = 344; % J*m^-2
duP = 2.5e-3;
maxCMOD = 0.1; % mm

[CMOD, P, u, nodeEla, elemEla] = SEBFunctionV2('bilinear', sigma_c, G_c, duP, maxCMOD, elems, nodes, 0.001);

%% Stress
dispV = reshape(u(1 : size(nodeEla, 1)*2, end), 2, [])';
% uy
figure
axis equal;
PlotContour(nodeEla,elemEla,dispV(:, 2),'uy',1);
axis off;

% ux
figure
axis equal;
PlotContour(nodeEla,elemEla,dispV(:, 1),'ux',1);
axis off;

vtkName = ['test','.vtk'];
vtkwrite(vtkName, 'UNSTRUCTURED_GRID', nodeEla(:,1), nodeEla(:,2), nodeEla(:,1)*0, ...
    'cells', elemEla, 'cell_types', 9, 'vectors', 'u', dispV(:, 1)', dispV(:, 2)', 0*dispV(:, 1)');
