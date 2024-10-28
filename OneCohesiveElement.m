% one cohesive element

clear; close all; clc
addpath("Func\")
%%  ***  Reas Ansys Mesh  ***

node = [0,0; 2,0; 2,0; 0,0];
nodeBou = {[1,0,0; 2,2,0]; [2,2,0; 3,2,0]; [3,2,0; 4,0,0]; [4,0,0; 1,0,0]};
elem = [1,2,3,4];
sumNode = size(node,1);

%% ***  Material para  *** (Ambati's Paper)
E = 0;

u_0 = zeros(numel(node), 1);

% Set boundary condition.
ubar = -1e-5;
fixNode = [];
nodeForce = [];

% (0, 0)
m = 1;
[fixNode] = [fixNode; m, 1, 0; m, 2, 0];

% (2, 0)
m = 2;
[fixNode] = [fixNode; m, 1, 0; m, 2, 0];

% line3 Y
line = nodeBou{3};
[fixNode] = [fixNode; line(:, 1), 2*ones(size(line, 1), 1), ubar * ones(size(line, 1), 1)];

fixNode = sortrows(fixNode, 1);

BC = setBC(fixNode, nodeForce, sumNode*2);

for i = 1 : size(BC.DirchletDOF, 1)
    ind = BC.DirchletDOF(i);
    u_0(ind) = BC.Dirichlet(i);
end

Para.ndim = 2; % dim
Para.isStress = 2;  % 1 - plane stress, 2 - plane strain
Para.E = E; % Young's Modulus based on (N/mm2)
Para.nu = 0.35; % Poisson's Ratio
Para.lambda = Para.E*Para.nu/((1+Para.nu)*(1-2*Para.nu)); % Lame Constant
Para.mu = Para.E/(2*(1+Para.nu)); % Lame Constant
Para.NNd = size(node,1); % number of nodes

node = node(:, 1 : Para.ndim);
u_c0 = reshape(u_0, 2, [])';
[GaussInfo] = shapeFunc_valueDeriv_CZM(elem, node, u_c0, Para);

%% Elastic problem with Cohesive zone model
sigma_c = 3.56; % MPa
G_c = 344; % J*m^-2
delta_c = 2*G_c/sigma_c*1e-3; % mm
lamda_cr = 0.001;

% Newton-Raphson method
KCZM_0 = globalK2D_CZM_Bilinear(Para, elem, GaussInfo, u_c0, 'bilinear', sigma_c, lamda_cr, delta_c);
Fc_0 = nodeForce_CZM(Para, elem, GaussInfo, u_c0, 'bilinear', sigma_c, lamda_cr, delta_c);

R_0 = Fc_0 - BC.RHS;

% Assemble the total K.
dRdu_0 = KCZM_0;

indIter = 0;
tole = 1e-5;
r = max(abs(R_0));
while max(abs(R_0)) > tole && indIter<100
    % solve
    % Big number
    max_abs_value = max(max(abs(dRdu_0)));
    n = floor(log10(max_abs_value));
    bigN = 10^(7+n);
    for i = 1 : size(BC.DirchletDOF, 1)
        ind = BC.DirchletDOF(i);
        R_0(ind) = 0 * dRdu_0(ind, ind) * bigN;
        dRdu_0(ind, ind) = dRdu_0(ind, ind) * bigN;
    end
    du_0 = -dRdu_0\R_0;

    u_1 = u_0 + du_0;
    u_c1 = reshape(u_1, 2, [])';
    [GaussInfo] = shapeFunc_valueDeriv_CZM(elem, node, u_c1, Para);
    KCZM_0 = globalK2D_CZM_Bilinear(Para, elem, GaussInfo, u_c1, 'bilinear', sigma_c, lamda_cr, delta_c);
    Fc_0 = nodeForce_CZM(Para, elem, GaussInfo, u_c1, 'bilinear', sigma_c, lamda_cr, delta_c);
    R_0 = - BC.RHS;
    R_0 = Fc_0 + R_0;
    u_0 = u_1;

    % Assemble the total K.
    dRdu_0 = KCZM_0;

    indIter = indIter + 1;
    [r] = [r; max(abs(R_0))];
end


%% Plot
% uy
figure
axis equal;
PlotContour(node{1},elem{1},du_0(1:size(node{1}, 1), 2),'uy',1);
axis off;

% ux
figure
axis equal;
PlotContour(node{1},elem{1},du_0(1:size(node{1}, 1), 1),'ux',1);
axis off;