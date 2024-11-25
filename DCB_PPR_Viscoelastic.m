% DCB with PPR CZM
% 基于粘弹性内聚模型的固体发动机界面力学性能研究_崔辉如.pdf

clear; clc; close;
addpath("Func\")

a = 150;
b = 5;
O = [0,0;];
numX = 150;
numY = 5;
% numX = 3;
% numY = 2;
isCoh = 0;
[node, elem] = generateMeshFEM('rect', a, b, O, numX, numY, isCoh);

% Find cohesive node
coNodesUp = node((node(:, 2) >= 20) & (node(:, 3) == 0), :);
coNodesUp = sortrows(coNodesUp, 2);
coNodesDown = [(1 : size(coNodesUp, 1))' + size(node, 1), coNodesUp(:, 2 : 3)];

% Cohesive element
coElme = [];
for i = 1 : size(coNodesUp, 1) - 1
    [coElme] = [coElme; coNodesUp(i, 1), coNodesUp(i + 1, 1), ...
        coNodesDown(i + 1, 1), coNodesDown(i, 1)];
end
[elem] = [elem; (1 : size(coElme, 1))' + size(elem, 1), 2 * ones(size(coElme, 1), 1), coElme];

[node] = [node; coNodesDown];

elemCoh = elem(elem(:, 2) == 2, 3 : end);
elemEla = elem(elem(:, 2) == 1, 3 : end);
node = node(:, 2 : 3);
elem = elem(:, 3 : end);

%% PPR cohesive model parameters
sigma_max_0 = 0.5; % MPa
alpha_k = 3.10;
delta_n = 2.5; % mm
tau_max = 0.5; % MPa
beta_k = 3.10;
delta_t = 2.5; % mm
m0 = 0.1672;
n0 = m0;

figure

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
delatT = 0.1;
tMax = 10;
vec = 0.48; % mm/s
[DMax] = initial_DMax(elem, Para);
P = zeros(tMax / delatT + 1, 1);
CMOD = zeros(tMax / delatT + 1, 1);
q = 1;
while t <= tMax
    % Boundary condition
    uS = vec * t;

    % A starting value for the unknown must be chosen;
    % usually the solution u_n from the last time step (n) is selected
    % modify u_n to satisfy the boundary conditions.

    % Step of the iteration: v = 0
    [BC, uv] = startingValue(un, uS, node);

    % Calculate the cohesive nodal forces and the stiffness matrix
    [KCZM, Fcv] = CohesiveMatrix(uv, elemCoh, node, Para, ParaPPR, DMax);

    % Residual
    Rv = K * uv + Fcv - BC.RHS;
    % Tanget matrix
    dRduv = K + KCZM;

    tempP = - (Rv(1) + Rv(3));

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

        tempP =  - (Rv(1) + Rv(3));

        % Modify Rv and dRduv to satisfy the boundary conditions.
        [Rv, dRduv] = ApplyBC(BC, Rv, dRduv);

        idxIter = idxIter + 1;
    end

    P(q) = tempP;
    CMOD(q) = uS;
    fprintf('Time %d, iteration number %d, uP = %f, P = %f\n', t, idxIter, uS, tempP);

    un = uv;
    t = t + delatT;
    q = q + 1;

    [~, ~, GaussInfo] = CohesiveMatrix(un, elemCoh, node, Para, ParaPPR, DMax);
    DMax = renew_DMax(Para, elemCoh, GaussInfo, reshape(un, 2, [])', DMax);

end
plot(CMOD, P / 2, 'LineWidth', 1.5)
xlim([-1 3])
ylim([-0.3 0.3])
grid on
hold on

xtitle = 'Separation/mm';
ytitle = 'Reaction force/N';
lgd = legend('m = 0.1677', 'm = 0.3344', 'm = 0.5066');
setPlotV2(xtitle, ytitle, lgd)

%% Sub function
function [BC, uv] = startingValue(un, uP, node)
uv = un;

% Boundary condition
fixNode = [];
nodeForce = [];

% (-1, 0) and (1, 0)
[fixNode] = [fixNode; 1, 1, 0];
% [fixNode] = [fixNode; 1, 2, 0];
[fixNode] = [fixNode; 2, 1, 0; 2, 2, 0];

% (0, 5) P
[fixNode] = [fixNode; 3, 1, uP; 4, 1, uP];

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