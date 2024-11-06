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
%%  ***  Read Abaqus Mesh  ***
parts = loadinp(YourModel);

% Copy the half model
nodes = [(1 : size(parts.Nodes, 1))'; parts.Nodes];

elmes = [];
for i = 1 : size(parts.Elements, 1)
    [elmes] = [elmes; i, 1, parts.Elements(i).Connectivity];
end

insertionRegions = [0, 0, 0, 19];
[newNodes, newElements, cohesiveElements] = insertCohesiveElements(nodes, elmes, insertionRegions);
