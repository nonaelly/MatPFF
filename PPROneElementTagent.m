% One element with PPR CZM

clear; clc; close;
node = [1,-1,0; 2,1,0; 3,1,0; 4,-1,0];
elem = [1, 2, 1, 2, 3, 4];
elemCoh = elem(elem(:, 2) == 2, 3 : end);
node = node(:, 2 : 3);
elem = elem(:, 3 : end);
%% PPR cohesive model parameters
% sigma_max = 0.5; % MPa
sigma_max = 0.25; % MPa
mA = [0.1677, 0.3344, 0.5066];
alpha_k = 3.10;
delta_n = 2.5; % mm
% tau_max = 0.5; % MPa
tau_max = 0.25; % MPa
nA = mA;
beta_k = 3.10;
delta_t = 2.5; % mm

figure
for k = 1 : 3
    m = mA(k);
    n = nA(k);
    ParaPPR = [delta_n, delta_t, m, n, alpha_k, beta_k, sigma_max, tau_max];
    %% Material para
    E = 1;
    Para.ndim = 2; % dim
    Para.isStress = 2;  % 1 - plane stress, 2 - plane strain
    Para.E = E; % Young's Modulus based on (N/mm2)
    Para.nu = 0.3; % Poisson's Ratio
    Para.lambda = Para.E*Para.nu/((1+Para.nu)*(1-2*Para.nu)); % Lame Constant
    Para.mu = Para.E/(2*(1+Para.nu)); % Lame Constant
    Para.NNd = size(node,1); % number of nodes

    %% Newton-Raphson method

    tole = 1e-6;

    % Separation
    uS = 0;
    un = zeros(2*Para.NNd, 1);
    t = 0;
    [DMax] = initial_DMax(elem, Para);
    P = zeros(301, 1);
    CMOD = zeros(301, 1);
    q = 1;
    while t <= 300
        % Boundary condition
        uS = - 1 / 100 * t * (t<=100) + (1 / 50 * (t - 100) - 1) * (t>100 && t<=150)...
            + (1 / 50 * (t - 150)) * (t>150 && t<=200) + (1.5 / 100 * (t - 200) + 1) * (t>200 && t<=300);

        % A starting value for the unknown must be chosen;
        % usually the solution u_n from the last time step (n) is selected
        % modify u_n to satisfy the boundary conditions.

        % Step of the iteration: v = 0
        [BC, uv] = startingValue(un, uS, node);

        % Calculate the cohesive nodal forces and the stiffness matrix
        [KCZM, Fcv] = CohesiveMatrix(uv, elemCoh, node, Para, ParaPPR, DMax);

        % Residual
        Rv = Fcv - BC.RHS;
        % Tanget matrix
        dRduv = KCZM;

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
            Rv = Fcv - BC.RHS;
            % Tanget matrix
            dRduv = KCZM;

            tempP =  - (Rv(1) + Rv(3));

            % Modify Rv and dRduv to satisfy the boundary conditions.
            [Rv, dRduv] = ApplyBC(BC, Rv, dRduv);

            idxIter = idxIter + 1;
        end

        P(q) = tempP;
        CMOD(q) = uS;
        fprintf('Time %d, iteration number %d, uP = %f, P = %f\n', t, idxIter, uS, tempP);

        un = uv;
        t = t + 0.2;
        q = q + 1;

        [~, ~, GaussInfo] = CohesiveMatrix(un, elemCoh, node, Para, ParaPPR, DMax);
        DMax = renew_DMax(Para, elemCoh, GaussInfo, reshape(un, 2, [])', DMax);

    end
    plot(CMOD, P / 2, 'LineWidth', 1.5)
    xlim([-1 3])
    ylim([-0.3 0.3])
    grid on
    hold on
end

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

% (-1, 0) and (1, 0) P
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