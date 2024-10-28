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

dx = 2;
numY = [19, 81/1.5; 0, 81/1.5];
% dx = 1;
% numY = [19, 81/1; 0, 81/1];
% dx = 10;
% numY = [10, 3; 0, 3];
[node, elem, nodeBou, elemBou] = deal(cell(2, 1));
isC = [0, 1];
sumNode = zeros(2, 1);
for i = 1:2
    [node{i}, elem{i}, nodeBou{i}, elemBou{i}] = generateMeshFEM('SEB', dx, numY(i, :), isC(i));
    sumNode(i) = size(node{i},1);
    node{i}(:,1)   = [];
    elem{i}(:,1:2) = [];
end

%% ***  Material para  *** (Ambati's Paper)
E = [14.2e3, 0];
[Para, GaussInfo] = deal(cell(2, 1));

u_0 = zeros(numel(node{1}) + numel(node{2})/2, 1);

% Set boundary condition.
ubar = -1e-5;
[fixNode, nodeForce] = generateBC_CZM_Bilinear(nodeBou, node, ubar);
BC = setBC(fixNode, nodeForce, (sumNode(1) + sumNode(2)/2)*2);

for i = 1 : size(BC.DirchletDOF, 1)
    ind = BC.DirchletDOF(i);
    u_0(ind) = BC.Dirichlet(i);
end

% (0, 19)
m = find(ismember(node{1}, [0, 19], 'rows'));
u_c0 = reshape(u_0([1+numel(node{1}):end, 2*m-1 : 2*m-2 + numel(node{2})/2]), 2, [])';

for i = 1:2
    Para{i}.ndim = 2; % dim
    Para{i}.isStress = 2;  % 1 - plane stress, 2 - plane strain
    Para{i}.E = E(i); % Young's Modulus based on (N/mm2)
    Para{i}.nu = 0.35; % Poisson's Ratio
    Para{i}.lambda = Para{i}.E*Para{i}.nu/((1+Para{i}.nu)*(1-2*Para{i}.nu)); % Lame Constant
    Para{i}.mu = Para{i}.E/(2*(1+Para{i}.nu)); % Lame Constant
    Para{i}.NNd = size(node{i},1); % number of nodes

    node{i} = node{i}(:, 1 : Para{i}.ndim);

    if i == 2
        [GaussInfo{i}] = shapeFunc_valueDeriv_CZM(elem{i}, node{i}, u_c0, Para{i});
    else
        [GaussInfo{i}] = shapeFunc_valueDeriv(elem{i}, node{i}, Para{i});
    end
end

K = globalK2D(Para{1}, elem{1}, GaussInfo{1});

%% Elastic problem with Cohesive zone model
sigma_c = 3.56; % MPa
G_c = 344; % J*m^-2
delta_c = 2*G_c/sigma_c*1e-3; % mm
lamda_cr = 0.001;

% Newton-Raphson method
KCZM_0 = globalK2D_CZM_Bilinear(Para{2}, elem{2}, GaussInfo{2}, u_c0, 'bilinear', sigma_c, lamda_cr, delta_c);
Fc_0 = nodeForce_CZM(Para{2}, elem{2}, GaussInfo{2}, u_c0, 'bilinear', sigma_c, lamda_cr, delta_c);
R_0 = - BC.RHS;
R_0(1:size(K, 1)) = K * u_0(1:size(K, 1)) + R_0(1:size(K, 1));
R_0([1+size(K, 1) : end, 2*m-1 : 2*m-2 + size(KCZM_0, 1)/2]) = Fc_0 + R_0([1+size(K, 1) : end, 2*m-1 : 2*m-2 + size(KCZM_0, 1)/2]);

% Assemble the total K.
dRdu_0 = zeros(size(K, 1) + size(KCZM_0, 1)/2);
idx =  [size(K, 1) + 1 : size(K, 1) + size(KCZM_0, 1)/2, 2*m-1 : 2*m-2 + size(KCZM_0, 1)/2];
dRdu_0(idx, idx) = dRdu_0(idx, idx) + KCZM_0;
idx =  1 : size(K, 1);
dRdu_0(idx, idx) = dRdu_0(idx, idx) + K;

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
    u_c1 = reshape(u_1([1+numel(node{1}):end, 2*m-1 : 2*m-2 + numel(node{2})/2]), 2, [])';
    [GaussInfo{2}] = shapeFunc_valueDeriv_CZM(elem{2}, node{2}, u_c1, Para{2});
    KCZM_0 = globalK2D_CZM_Bilinear(Para{2}, elem{2}, GaussInfo{2}, u_c1, 'bilinear', sigma_c, lamda_cr, delta_c);
    Fc_0 = nodeForce_CZM(Para{2}, elem{2}, GaussInfo{2}, u_c1, 'bilinear', sigma_c, lamda_cr, delta_c);
    R_0 = - BC.RHS;
    R_0(1:size(K, 1)) = K * u_1(1:size(K, 1)) + R_0(1:size(K, 1));
    R_0([1+size(K, 1) : end, 2*m-1 : 2*m-2 + size(KCZM_0, 1)/2]) = Fc_0 + R_0([1+size(K, 1) : end, 2*m-1 : 2*m-2 + size(KCZM_0, 1)/2]);
    u_0 = u_1;

    % Assemble the total K.
    dRdu_0 = zeros(size(K, 1) + size(KCZM_0, 1)/2);
    idx =  [size(K, 1) + 1 : size(K, 1) + size(KCZM_0, 1)/2, 2*m-1 : 2*m-2 + size(KCZM_0, 1)/2];
    dRdu_0(idx, idx) = dRdu_0(idx, idx) + KCZM_0;
    idx =  1 : size(K, 1);
    dRdu_0(idx, idx) = dRdu_0(idx, idx) + K;

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