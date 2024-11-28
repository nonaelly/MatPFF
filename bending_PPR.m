% mixed-mode bending test with PPR CZM
% 'Computational implementation of the PPR potential-based cohesive model
%   in ABAQUS: Educational perspective'

clear; clc; close all;
addpath("Func\")

L = 51;
h = 1.56;
O = [0, 0; 0, -h];
numX = 2 * L / 0.1;
numY = 10;
% numX = 3;
% numY = 2 * 2;
isCoh = 0;
node = [];
elem = [];
[nodePart, elemPart] = deal(cell(2, 1));
for i = 1 : 2
    [nodePart{i}, elemPart{i}] = generateMeshFEM('rect', 2 * L, h, O(i, :), numX, numY, isCoh);
    if i > 1
        nodePart{i}(:, 1) = nodePart{i - 1}(end, 1) + nodePart{i}(:, 1);
        elemPart{i}(:, 1) = elemPart{i - 1}(end, 1) + elemPart{i}(:, 1);
        elemPart{i}(:, 3 : end) = nodePart{i - 1}(end, 1) + elemPart{i}(:, 3 : end);
    end
    [node] = [node; nodePart{i}];
    [elem] = [elem; elemPart{i}];
end

a0 = 33.7;
% Find cohesive node
coNodesUp = nodePart{1}((nodePart{1}(:, 2) <= 2 * L - a0) & (nodePart{1}(:, 3) == 0), :);
coNodesUp = sortrows(coNodesUp, 2);

coNodesDown = nodePart{2}((nodePart{2}(:, 2) <= 2 * L - a0) & (nodePart{2}(:, 3) == 0), :);
coNodesDown = sortrows(coNodesDown, 2);

% Cohesive element
coElme = [];
for i = 1 : size(coNodesUp, 1) - 1
    [coElme] = [coElme; coNodesDown(i, 1), coNodesDown(i + 1, 1), ...
        coNodesUp(i + 1, 1), coNodesUp(i, 1)];
end
[elem] = [elem; (1 : size(coElme, 1))' + size(elem, 1), 2 * ones(size(coElme, 1), 1), coElme];

elemCoh = elem(elem(:, 2) == 2, 3 : end);
elemEla = elem(elem(:, 2) == 1, 3 : end);
nodeEla = node(1 : max(max(elemEla)), 2 : 3);

node = node(:, 2 : 3);
elem = elem(:, 3 : end);

%% PPR cohesive model parameters
phi_n = 500;
lambda_n = 0.02;
sigma_max = 100; % MPa
alpha_k = 3;


m = -(lambda_n^2 * alpha_k * (alpha_k - 1)) / (lambda_n^2 * alpha_k - 1);
delta_n = phi_n * (alpha_k * lambda_n * (1 - lambda_n)^(alpha_k - 1) * ...
    (alpha_k / m + 1) * (alpha_k / m * lambda_n + 1)^(m-1)) / sigma_max * 1e-3; % mm

tau_max = sigma_max; % MPa
beta_k = alpha_k;
delta_t = delta_n; % mm
n = m;

ParaPPR = [delta_n, delta_t, m, n, alpha_k, beta_k, sigma_max, tau_max];
%% Material para
Para.ndim = 2; % dim
Para.isStress = 1;  % 1 - plane stress, 2 - plane strain
Para.E = 122e3; % Young's Modulus based on (N/mm2)
Para.nu = 0.25; % Poisson's Ratio
Para.lambda = Para.E*Para.nu/((1+Para.nu)*(1-2*Para.nu)); % Lame Constant
Para.mu = Para.E/(2*(1+Para.nu)); % Lame Constant
Para.NNd = size(node,1); % number of nodes

%% Generate elastic stiffness matrix
[GaussInfoEla] = shapeFunc_valueDeriv(elemEla, node, Para);
K = globalK2D(Para, elemEla, GaussInfoEla);

%% Newton-Raphson method

tole = 1e-6;

% Separation
PLever = 0;
un = zeros(2*Para.NNd, 1);
t = 0;
delatT = 0.02;
vec = 1; % mm/s
tMax = 5/vec;

[DMax] = initial_DMax(elem, Para);
P = [];
CMOD = [];
q = 1;
u = [];
numVtk = 0;
rangeVtk = 1;
idxD = find(ismember(node, [2 * L, 0], 'rows'));

while t <= tMax
    % Boundary condition
    PLever = vec * t;

    % A starting value for the unknown must be chosen;
    % usually the solution u_n from the last time step (n) is selected
    % modify u_n to satisfy the boundary conditions.

    % Step of the iteration: v = 0
    [BC, uv] = startingValue(un, PLever, node, coNodesDown);

    % Calculate the cohesive nodal forces and the stiffness matrix
    [KCZM, Fcv] = CohesiveMatrix(uv, elemCoh, node, Para, ParaPPR, DMax);

    % Residual
    Rv = K * uv + Fcv - BC.RHS;
    % Tanget matrix
    dRduv = K + KCZM;

    % Modify Rv and dRduv to satisfy the boundary conditions.
    [Rv, dRduv] = ApplyBC(BC, Rv, dRduv);

    idxIter = 0;

    while max(abs(Rv)) > tole && idxIter < 20
        % step of the iteration: v = 1, 2, ...
        duv = -dRduv\Rv;

        % step of the iteration: v + 1, renew u, R, dRdu
        uv = uv + duv;

        % Calculate the cohesive nodal forces and the stiffness matrix
        [KCZM, Fcv] = CohesiveMatrix(uv, elemCoh, node, Para, ParaPPR, DMax);

        % Residual
        Rv = K * uv + Fcv - BC.RHS;
        % Tanget matrix
        dRduv = K + KCZM;

        % Modify Rv and dRduv to satisfy the boundary conditions.
        [Rv, dRduv] = ApplyBC(BC, Rv, dRduv);

        idxIter = idxIter + 1;
    end

    un = uv;
    dispV = reshape(un(1 : size(nodeEla, 1)*2), 2, [])';
    t = t + delatT;
    [u] = [u, un];
    q = q + 1;

    [~, ~, GaussInfo] = CohesiveMatrix(un, elemCoh, node, Para, ParaPPR, DMax);
    DMax = renew_DMax(Para, elemCoh, GaussInfo, reshape(un, 2, [])', DMax);

    DLever = un(idxD(1) * 2);
    [P] = [P, PLever];
    [CMOD] = [CMOD, DLever];
    fprintf('Time %f, iteration number %d, d = %f, P = %f\n', t, idxIter, DLever, PLever);

    if idxIter == 20
        1;
    end
    %     if mod(t, rangeVtk) < 1e-8 || mod(t, rangeVtk) > rangeVtk - 1e-8
    %         dispV = reshape(un(1 : size(nodeEla, 1)*2), 2, [])';
    %         numVtk = numVtk + 1;
    %         vtkName = ['PPR-test-complete-model-', num2str(numVtk),'.vtk'];
    %         vtkwrite(vtkName, 'UNSTRUCTURED_GRID', nodeEla(:,1), nodeEla(:,2), nodeEla(:,1)*0, ...
    %             'cells', elemEla, 'cell_types', 9, 'vectors', 'u', ...
    %             dispV(:, 1)', dispV(:, 2)', 0*dispV(:, 1)');
    %     end
end
%%
figure
plot(CMOD, P, 'LineWidth', 1.5)
lgd = legend('line1');

% xlim([-1 3])
% ylim([-0.3 0.3])
grid on
hold on

xtitle = 'Separation/mm';
ytitle = 'Reaction force/N';
setPlotV2(xtitle, ytitle, lgd)

%% Vtk

dispV = reshape(un(1 : size(nodeEla, 1)*2), 2, [])';
numVtk = numVtk + 1;
vtkName = ['PPR-bending-model-', num2str(numVtk),'.vtk'];
vtkwrite(vtkName, 'UNSTRUCTURED_GRID', nodeEla(:,1), nodeEla(:,2), nodeEla(:,1)*0, ...
    'cells', elemEla, 'cell_types', 9, 'vectors', 'u', ...
    dispV(:, 1)', dispV(:, 2)', 0*dispV(:, 1)');

%% Sub function
function [BC, uv] = startingValue(un, P, node, coNodesDown)
uv = un;

% Boundary condition
fixNode = [];
nodeForce = [];

c = 60;
L = 51;

% (L, h)
idxP = find(ismember(node, [51, 1.56], 'rows'));
[nodeForce] = [nodeForce; idxP, 2, - P * (c + L) / L];

% (2L, h)
idxP = find(ismember(node, [51 * 2, 1.56], 'rows'));
[nodeForce] = [nodeForce; idxP, 2, P * c / L];

% (0, -h) uy = 0
idxF = find(ismember(node, [0, -1.56], 'rows'));
[fixNode] = [fixNode; idxF, 2, 0];

% (2L, -h) ux = uy = 0
idxF = find(ismember(node, [2 * 51, -1.56], 'rows'));
[fixNode] = [fixNode; idxF, 2, 0; idxF, 1, 0];

fixNode = sortrows(fixNode, 1);

BC = setBC(fixNode, nodeForce, size(node, 1) * 2);

% uv should satisfy the boundary condition
uv(BC.DirchletDOF) = BC.Dirichlet;
end

function [KCZM, Fcv, GaussInfo] = CohesiveMatrix(uv, elemCoh, node, Para, ParaPPR, DMax)

uv = reshape(uv, 2, [])';
GaussInfo = shapeFunc_valueDeriv_CZM(elemCoh, node, uv, Para);
[KCZM, Fcv] = globalK2D_CZM_PPR(Para, elemCoh, GaussInfo, uv, 'PPR', ParaPPR, DMax);
end

function [Rv, dRduv] = ApplyBC(BC, Rv, dRduv)
% Because Δu = - dRduv \ Rv
Rv(BC.DirchletDOF) = 0;
dRduv(:, BC.DirchletDOF) = 0;
dRduv(BC.DirchletDOF, :) = 0;
for i = 1 : size(BC.DirchletDOF)
    dRduv(BC.DirchletDOF(i), BC.DirchletDOF(i)) = 1;
end
end