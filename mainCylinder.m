% Elastic problem.
% Thick cylinder

clear; close all; clc
addpath("Func\")
%%  ***  Reas Ansys Mesh  ***
R1 = 5;
R2 = 10;
O = [0,0];
numR = 12;
numTheta = 12;
% numR = 1;
% numTheta = 1;
[node, elem, nodeBou, elemBou] = generateMeshFEM('cyl', R1, R2, O, numR, numTheta);
sumNode = size(node,1);

%% ***  Material para  *** (Ambati's Paper)
% Para.PFModel = 2; % 1-AT2; 2-AT1
Para.ndim = 2; % dim
Para.isStress = 2;  % 1 - plane stress, 2 - plane strain
Para.E = 25714; % Young's Modulus based on (N/mm2)
Para.nu = 0.2857; % Poisson's Ratio
Para.lambda = Para.E*Para.nu/((1+Para.nu)*(1-2*Para.nu)); % Lame Constant
Para.mu = Para.E/(2*(1+Para.nu)); % Lame Constant
% Para.Gc = 0.089; % Critical energy release for unstable crack (Gc, N/mm)
% Para.Len = 3; %

Para.NNd = size(node,1); % number of nodes

elem(:,1:2) = [];
node(:,1)   = [];
node = node(:, 1 : Para.ndim);

[GaussInfo] = shapeFunc_valueDeriv(elem, node, Para);

%% Elastic problem

K = globalK2D(Para, elem, GaussInfo);

% Set boundary condition.
boundary = {'p', 'dy', 'f', 'dx'};
pressure = 1;
[fixNode, nodeForce] = generateBC(boundary, nodeBou, elemBou, elem, node, pressure);

BC = setBC(fixNode, nodeForce, sumNode*2);

% solve
F = BC.RHS;
% Big number
max_abs_value = max(max(abs(K)));
n = floor(log10(max_abs_value));
bigN = 10^(9+n);
for i = 1 : size(BC.DirchletDOF, 1)
    ind = BC.DirchletDOF(i);
    F(ind) = BC.Dirichlet(i) * K(ind, ind) * bigN;
    K(ind, ind) = K(ind, ind) * bigN;
end
Disp = K\F;
Stress = calcStress2DV2(GaussInfo, elem, Para, Disp);

Disp = reshape(Disp,Para.ndim,[]);
Disp = Disp';

% ux
figure
axis equal;
PlotContour(node,elem,Disp(:,1),'ux',1);
axis off;
% uy
figure
axis equal;
PlotContour(node,elem,Disp(:,2),'uy',1);
axis off;
% mises
figure
axis equal;
PlotContour(node,elem,Stress.vonMises,'mises',1);
axis off;
