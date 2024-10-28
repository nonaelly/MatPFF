% single-edge notched beam (SE(B)) test
% <A bilinear cohesive zone model tailored for fracture of asphalt concrete
%   considering viscoelastic bulk material>

clear; close all; clc
addpath("Func\")
%%  ***  Reas Ansys Mesh  ***
w = 100;
a = 19;
t = 75;
b1 = 162;
b2 = 26;

% dx = 2;
% numY = [19, 81/1.5; 0, 81/1.5];
dx = 1;
numY = [19, 81/1; 0, 81/1];
[node, elem, nodeBou, elemBou] = deal(cell(2, 1));
isC = [0, 1];
sumNode = zeros(2, 1);
for i = 1:2
    [node{i}, elem{i}, nodeBou{i}, elemBou{i}] = generateMeshFEM('SEB', dx, numY(i, :), isC(i));
    sumNode(i) = size(node{i},1);
end

%% ***  Material para  *** (Ambati's Paper)
E = [14.2e3, 0];
[Para, GaussInfo] = deal(cell(2, 1));
for i = 1:2
    Para{i}.ndim = 2; % dim
    Para{i}.isStress = 2;  % 1 - plane stress, 2 - plane strain
    Para{i}.E = E(i); % Young's Modulus based on (N/mm2)
    Para{i}.nu = 0.35; % Poisson's Ratio
    Para{i}.lambda = Para{i}.E*Para{i}.nu/((1+Para{i}.nu)*(1-2*Para{i}.nu)); % Lame Constant
    Para{i}.mu = Para{i}.E/(2*(1+Para{i}.nu)); % Lame Constant
    Para{i}.NNd = size(node{i},1); % number of nodes

    elem{i}(:,1:2) = [];
    node{i}(:,1)   = [];
    node{i} = node{i}(:, 1 : Para{i}.ndim);
    if i == 2
        u_0 = zeros(size(node{i}));
        [GaussInfo{i}] = shapeFunc_valueDeriv_CZM(elem{i}, node{i}, u_0, Para{i});
    else
        [GaussInfo{i}] = shapeFunc_valueDeriv(elem{i}, node{i}, Para{i});
    end
end

K = globalK2D(Para{1}, elem{1}, GaussInfo{1});

%% Elastic problem with Cohesive zone model
sigma_c = 3.56; % MPa
G_c = 344; % J*m^-2
delta_c = 2*G_c/sigma_c*1e-3; % mm
lamda_cr = 0.04;

% Newton-Raphson method
Tc_0 = zeros(numel(node{2}), 1);
Fc_0 = nodeForce_CZM(Tc_0);
R_0 = K * u_0 + Fc_0 - BC.RHS;
KCZM_0 = globalK2D_CZM_Bilinear(Para{2}, elem{2}, GaussInfo{2}, u_0, 'bilinear', sigma_c, lamda_cr, delta_c);
dRdu_0 = KCZM_0 + K;
indIter = 0;
while R_0 > tole && indIter<100
du_0 = -dRdu_0\R_0;
u_1 = u_0 + du_0;
dTc_0 = Dc * Bc * du_0;
Tc_1 = dTc_0 + Tc_0;
Fc_1 = nodeForce_CZM(Tc_1);
R_1 = K * u_1 + Fc_1 - BC.RHS;
KCZM_1 = globalK2D_CZM_Bilinear(Para{2}, elem{2}, GaussInfo{2}, u_1, 'bilinear', sigma_c, lamda_cr, delta_c);
dRdu_1 = KCZM_1 + K;
du_1 = -dRdu_1\R_1;
indIter = indIter + 1;
end

% Assemble the total K.
Ktotal = zeros(size(K, 1) + size(KCZM_0, 1)/2);

idx =  [size(K, 1) + 1 : size(K, 1) + size(KCZM_0, 1)/2, 1 : size(KCZM_0, 1)/2];
Ktotal(idx, idx) = Ktotal(idx, idx) + KCZM_0;
idx =  1 : size(K, 1);
Ktotal(idx, idx) = Ktotal(idx, idx) + K;

% Set boundary condition.
ubar = -1e-6;
[fixNode, nodeForce] = generateBC_CZM_Bilinear(nodeBou, node, ubar);

BC = setBC(fixNode, nodeForce, (sumNode(1) + sumNode(2)/2)*2);

% solve
Disp = solveBC(BC, Ktotal);
% Stress = calcStress2DV2(GaussInfo, elem, Para, Disp);

Disp = reshape(Disp,Para{1}.ndim,[]);
Disp = Disp';

%% Plot
% uy
figure
axis equal;
PlotContour(node{1},elem{1},Disp(1:size(node{1}, 1), 2),'uy',1);
axis off;

% ux
figure
axis equal;
PlotContour(node{1},elem{1},Disp(1:size(node{1}, 1), 1),'ux',1);
axis off;