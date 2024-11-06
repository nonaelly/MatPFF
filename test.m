clear; close all; clc
addpath("Func\")

%% Model 1
YourModel = 'Job-bulk2.inp';  % Choose your model
% YourModel = 'Job-CZM-mesh2.inp';  % Choose your model

sigma_c = 3.56; % MPa
G_c = 344; % J*m^-2
duP = 2.5e-3;
maxCMOD = 0.25; % mm

% lamda_cr = 0.001
%%  ***  Read Abaqus Mesh  ***
parts = loadinp(YourModel);

% Copy the half model
nodes = [(1 : size(parts.Nodes, 1))', parts.Nodes];

elmes = [];
for i = 1 : size(parts.Elements, 1)
    [elmes] = [elmes; i, 1, parts.Elements(i).Connectivity];
end

insertionRegions = [0, 0, 0, 19];
[newNodes, newElements, cohesiveElements] = insertCohesiveElements(nodes, elmes, insertionRegions);
elem = [newElements; cohesiveElements];
% Plot
figure(1);
hold on;
for i = 1:size(elem, 1)
    nodes = newNodes(elem(i, 3:end), 2:3);
    fill(nodes(:, 1), nodes(:, 2), 'w', 'EdgeColor', 'k');
    if elem(i, 2) == 2 % CZM
        plot(newNodes(elem(i, 3:end), 2), newNodes(elem(i, 3:end), 3), 'ro', 'MarkerFaceColor', 'r');
    end
end

title('single-edge notched beam (SE(B)) test');
xlabel('X');
ylabel('Y');
axis equal;

vtkName = ['test','.vtk'];
vtkwrite(vtkName, 'UNSTRUCTURED_GRID', newNodes(:,2), newNodes(:,3), newNodes(:,1)*0, ...
    'cells', elem(:, 3 : end), 'cell_types', 9);

