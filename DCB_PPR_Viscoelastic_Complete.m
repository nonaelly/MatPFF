% DCB with PPR CZM
% 基于粘弹性内聚模型的固体发动机界面力学性能研究_崔辉如.pdf

clear; clc; close all;
addpath("Func\")

a = 150;
b = 5;
O = [0, 0; 0, -5];
numX = 150;
numY = 5;
% numX = 3;
% numY = 2 * 2;
isCoh = 0;
node = [];
elem = [];
[nodePart, elemPart] = deal(cell(2, 1));
for i = 1 : 2
    [nodePart{i}, elemPart{i}] = generateMeshFEM('rect', a, b, O(i, :), numX, numY, isCoh);
    if i > 1
        nodePart{i}(:, 1) = nodePart{i - 1}(end, 1) + nodePart{i}(:, 1);
        elemPart{i}(:, 1) = elemPart{i - 1}(end, 1) + elemPart{i}(:, 1);
        elemPart{i}(:, 3 : end) = nodePart{i - 1}(end, 1) + elemPart{i}(:, 3 : end);
    end
    [node] = [node; nodePart{i}];
    [elem] = [elem; elemPart{i}];
end

% Find cohesive node
coNodesUp = nodePart{1}((nodePart{1}(:, 2) >= 20) & (nodePart{1}(:, 3) == 0), :);
coNodesUp = sortrows(coNodesUp, 2);

coNodesDown = nodePart{2}((nodePart{2}(:, 2) >= 20) & (nodePart{2}(:, 3) == 0), :);
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
elemIdxCoh = elem(elem(:, 2) == 2, 1);

node = node(:, 2 : 3);
elem = elem(:, 3 : end);

%% PPR cohesive model parameters
sigma_max_0 = 0.5; % MPa
alpha_k = 3.10;
delta_n = 2.5; % mm
m0 = 0.1672;

tau_max = sigma_max_0; % MPa
beta_k = alpha_k;
delta_t = delta_n; % mm
n0 = m0;

ParaPPR = [delta_n, delta_t, m0, n0, alpha_k, beta_k, sigma_max_0, tau_max];
%% Material para
Para.ndim = 2; % dim
Para.isStress = 2;  % 1 - plane stress, 2 - plane strain
Para.E = 72e3; % Young's Modulus based on (N/mm2)
Para.nu = 0.33; % Poisson's Ratio
Para.lambda = Para.E*Para.nu/((1+Para.nu)*(1-2*Para.nu)); % Lame Constant
Para.mu = Para.E/(2*(1+Para.nu)); % Lame Constant
Para.NNd = size(node,1); % number of nodes

%% Generate elastic stiffness matrix
[GaussInfoEla] = shapeFunc_valueDeriv(elemEla, node, Para);
K = globalK2D(Para, elemEla, GaussInfoEla);

%% Newton-Raphson method

tole = 1e-6;

% Separation
uS = 0;
un = zeros(2*Para.NNd, 1);
t = 0;
delatT = 0.03;
vec = 0.48; % mm/s
% tMax = 5 / vec;

% tMax = 5 / vec;
% 9.594
% tMax = 9.594;
tMax = 9.8;

[sVar] = initial_StateVaribles(elem, Para);
P = [];
CMOD = [];
q = 1;
u = [];
numVtk = 0;
rangeVtk = 1;
idxP = find(ismember(node, [0, 5], 'rows'));

while t <= tMax
    % Boundary condition
    uS = vec * t;

    % A starting value for the unknown must be chosen;
    % usually the solution u_n from the last time step (n) is selected
    % modify u_n to satisfy the boundary conditions.

    % Step of the iteration: v = 0
    [BC, uv] = startingValue(un, uS, node, coNodesDown);

    % Calculate the cohesive nodal forces and the stiffness matrix
    [KCZM, Fcv] = CohesiveMatrix(uv, elemCoh, node, Para, ParaPPR, sVar);

    % Residual
    Rv = K * uv + Fcv - BC.RHS;
    % Tanget matrix
    dRduv = K + KCZM;

    TempRv = K * uv + Fcv - BC.RHS;
    tempP = - TempRv(idxP * 2);

    % Modify Rv and dRduv to satisfy the boundary conditions.
    [Rv, dRduv] = ApplyBC(BC, Rv, dRduv);

    idxIter = 0;

    while max(abs(Rv)) > tole && idxIter < 20
        % step of the iteration: v = 1, 2, ...
        duv = -dRduv\Rv;

        % step of the iteration: v + 1, renew u, R, dRdu
        uv = uv + duv;

        % Calculate the cohesive nodal forces and the stiffness matrix
        [KCZM, Fcv] = CohesiveMatrix(uv, elemCoh, node, Para, ParaPPR, sVar);

        % Residual
        Rv = K * uv + Fcv - BC.RHS;
        % Tanget matrix
        dRduv = K + KCZM;

        % Modify Rv and dRduv to satisfy the boundary conditions.
        [Rv, dRduv] = ApplyBC(BC, Rv, dRduv);

        idxIter = idxIter + 1;

        if max(abs(Rv)) <= tole || idxIter >= 20
            tempFixP = K * uv;
            tempFixP = tempFixP(coNodesDown(end, 1) * 2);

            TempRv = K * uv + Fcv - BC.RHS;
            tempP = - TempRv(idxP * 2);
        end

    end

    [P] = [P, tempP];
    [CMOD] = [CMOD, uS];
    %     fprintf('Time %f, iteration number %d, uP = %f, P = %f\n', t, idxIter, uS, tempP);
    fprintf('Time %f, iteration number %d, uP = %f\n', t, idxIter, uS);

    un = uv;
    t = t + delatT;
    [u] = [u, un];
    q = q + 1;

    [~, ~, GaussInfo] = CohesiveMatrix(un, elemCoh, node, Para, ParaPPR, sVar);
    sVar = renew_DMax(Para, elem, elemIdxCoh, GaussInfo, reshape(un, 2, [])', sVar);

end

%% Plot
figure
plot(CMOD, -P, 'LineWidth', 1.5)
% nameLgd = cell(size(P, 1), 1);
% for i = 1 : size(P, 1)
%     plot(CMOD, -P(i, :), 'LineWidth', 1.5)
%     nameLgd{i} = ['line', num2str(i)];
%     hold on
% end
% lgd = legend(nameLgd);

% xlim([-1 3])
% ylim([-0.3 0.3])
grid on
hold on

xtitle = 'Separation/mm';
ytitle = 'Reaction force/N';
lgd = legend('v = 0.48 mm/s');
setPlotV2(xtitle, ytitle, lgd)
%% Vtk
dispV = reshape(un(1 : size(nodeEla, 1)*2), 2, [])';
numVtk = numVtk + 1;
vtkName = ['PPR-test-complete-model-', num2str(numVtk),'.vtk'];
vtkwrite(vtkName, 'UNSTRUCTURED_GRID', nodeEla(:,1), nodeEla(:,2), nodeEla(:,1)*0, ...
    'cells', elemEla, 'cell_types', 9, 'vectors', 'u', ...
    dispV(:, 1)', dispV(:, 2)', 0*dispV(:, 1)');

%% Stress
% % uy
% figure
% axis equal;
% PlotContour(nodeEla,elemEla,dispV(:, 2),'uy',1);
% axis off;
%
% % ux
% figure
% axis equal;
% PlotContour(nodeEla,elemEla,dispV(:, 1),'ux',1);
% axis off;


%% Sub function
function [BC, uv] = startingValue(un, uP, node, coNodesDown)
uv = un;

% Boundary condition
fixNode = [];
nodeForce = [];

% (150, 0)
idxF = find(ismember(node, [150, 0], 'rows'));
[fixNode] = [fixNode; idxF(1), 1, 0; idxF(2), 1, 0; idxF(1), 2, 0; idxF(2), 2, 0];
% [fixNode] = [fixNode; idxF(1), 1, 0; idxF(2), 1, 0];

% (0, 5) uP
idxP = find(ismember(node, [0, 5], 'rows'));
[fixNode] = [fixNode; idxP, 2, uP];

% (0, -5) -uP
idxP = find(ismember(node, [0, -5], 'rows'));
[fixNode] = [fixNode; idxP, 2, -uP];

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
    dRduv(BC.DirchletDOF(i), BC.DirchletDOF(i)) = 1e5;
end
end