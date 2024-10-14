% Elastic problem.
% Thick cylinder

clear; close all
addpath("Func\")
%%  ***  Reas Ansys Mesh  ***
R1 = 5;
R2 = 10;
O = [0,0];
% numR = 4;
% numTheta = 5;
numR = 1;
numTheta = 1;
[node, elem, nodeBou, elemBou] = generateMeshFEM('cyl', R1, R2, O, numR, numTheta);
sumNode = size(node,1);

%% ***  Material para  *** (Ambati's Paper)
% Para.PFModel = 2; % 1-AT2; 2-AT1
Para.ndim = 2; % dim
Para.isStress = 2;  % 1 - plane stress, 2 - plane strain
Para.E = 1; % Young's Modulus based on (N/mm2)
Para.nu = 0.3; % Poisson's Ratio
Para.lambda = Para.E*Para.nu/((1+Para.nu)*(1-2*Para.nu)); % Lame Constant
Para.mu = Para.E/(2*(1+Para.nu)); % Lame Constant
% Para.Gc = 0.089; % Critical energy release for unstable crack (Gc, N/mm)
% Para.Len = 3; %

Para.NNd = size(node,1); % number of nodes

elem(:,1:2) = [];
node(:,1)   = [];
node = node(:, 1 : Para.ndim);

[GaussInfo] = shapeFunc_valueDeriv(elem, node, Para);
% K = globalK2DV2(Para, elem, GaussInfo);

K = globalK2D(Para, elem, GaussInfo);

%% Elastic problem
% Set boundary condition.
boundary = {'p', 'dy', 'f', 'dx'};
fixNode = generateBC(boundary, nodeBou);

BC = ElasSENT(fixNode, sumNode*2, 0);

% solve
Disp = zeros(Para.NNd*Para.ndim,1);
F = zeros(Para.NNd*Para.ndim,1);
F(BC.FreeDOF) = F(BC.FreeDOF) - K(BC.FreeDOF, BC.DirchletDOF) * BC.Dirichlet;
Disp(BC.FreeDOF) = K(BC.FreeDOF, BC.FreeDOF) \ F(BC.FreeDOF);
Disp(BC.DirchletDOF) = BC.Dirichlet; % enforce BC 