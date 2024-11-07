% single-edge notched beam (SE(B)) test
% <A bilinear cohesive zone model tailored for fracture of asphalt concrete
%   considering viscoelastic bulk material>

clear; close all; clc
addpath("Func\")
%%  ***  Read Abaqus Mesh  ***

dx = 0.5;
numY = [19, 81/1.5];
% dx = 0.1;
% numY = [19, 81*2];
% dx = 0.1;
% numY = [19*2, 81*3];
% dx = 4;
% numY = [2, 3];
[node, elem, nodeBou, elemBou] = generateMeshFEM_SEB(dx, numY);

%% ***  Material para  *** (Ambati's Paper)
E = 14.2e3;
Para.ndim = 2; % dim
Para.isStress = 2;  % 1 - plane stress, 2 - plane strain
Para.E = E; % Young's Modulus based on (N/mm2)
Para.nu = 0.35; % Poisson's Ratio
Para.lambda = Para.E*Para.nu/((1+Para.nu)*(1-2*Para.nu)); % Lame Constant
Para.mu = Para.E/(2*(1+Para.nu)); % Lame Constant
Para.NNd = size(node,1); % number of nodes

%% Elastic problem with Cohesive zone model
sigma_c = 3.56; % MPa
G_c = 344; % J*m^-2
delta_c = 2*G_c/sigma_c*1e-3; % mm
lamda_cr = 0.001;

% Generate elastic stiffness matrix
GaussInfo = deal(cell(2, 1));
elemEla = elem(elem(:, 2) == 1, 3 : end);
[GaussInfo{1}] = shapeFunc_valueDeriv(elemEla, node(:, 2:3), Para);
K = globalK2D(Para, elemEla, GaussInfo{1});

% Newton-Raphson method
tole = 1e-6;

duP = 2.5e-3;
numStep = ceil(1/duP);
un = zeros(2*Para.NNd, numStep + 1);
CMOD = zeros(numStep, 1);
P = zeros(numStep, 1);
for n = 1 : numStep
    % Boundary condition
    uP = - n * duP;

    % A starting value for the unknown must be chosen;
    % usually the solution u_n from the last time step (n) is selected
    % modify u_n to satisfy the boundary conditions.

    % Step of the iteration: v = 0
    [BC, uv] = startingValue(un(:, n), uP, nodeBou, node);

    % Calculate the cohesive nodal forces and the stiffness matrix
    [KCZM, Fcv] = CohesiveMatrix(uv, sigma_c, lamda_cr, delta_c, elem, node, Para);

    % Residual
    Rv = K * uv + Fcv - BC.RHS;
    % Tanget matrix
    dRduv = K + KCZM;

    % Modify Rv and dRduv to satisfy the boundary conditions.
    [Rv, dRduv] = ApplyBC(BC, Rv, dRduv);

    idxIter = 0;

    while max(abs(Rv)) > tole && idxIter < 100
        % step of the iteration: v = 1, 2, ...
        duv = -dRduv\Rv;

        % step of the iteration: v + 1, renew u, R, dRdu
        uv = uv + duv;

        % Calculate the cohesive nodal forces and the stiffness matrix
        [KCZM, Fcv] = CohesiveMatrix(uv, sigma_c, lamda_cr, delta_c, elem, node, Para);

        % Residual
        Rv = K * uv + Fcv - BC.RHS;
        % Tanget matrix
        dRduv = K + KCZM;

        tempP = 2 * sum(Rv(BC.DirchletDOF(abs(BC.Dirichlet) > 0))) * 75 / 1e3;
%         tempP = (Rv(BC.DirchletDOF(1)));

        % Modify Rv and dRduv to satisfy the boundary conditions.
        [Rv, dRduv] = ApplyBC(BC, Rv, dRduv);

        idxIter = idxIter + 1;
    end

    un(:, n + 1) = uv;
    CMOD(n) = uv(1);
    P(n) = tempP;
    fprintf('Step %d, iteration number %d, uP = %f, P = %f, CMOD= %f\n', n, idxIter, uP, tempP, CMOD(n));
    if uv(1) > 0.25/2
        break
    end
%     if uv(1) > 0.01
%         break
%     end
end
figure
plot(2 * CMOD(1:n), -  P(1:n))
xlim([0 0.25])

%% Plot
elemEla = elem(elem(:, 2) == 1, 3 : end);
nodeEla = node(1 : max(max(elemEla)), 2:3);
dispV = reshape(uv(1 : size(nodeEla, 1)*2), 2, [])';
% uy
figure
axis equal;
PlotContour(nodeEla,elemEla,dispV(:, 2),'uy',1);
axis off;

% ux
figure
axis equal;
PlotContour(nodeEla,elemEla,dispV(:, 1),'ux',1);
axis off;

%% Sub function
function [BC, uv] = startingValue(un, uP, nodeBou, node)
uv = un;

% Boundary condition
[fixNode, nodeForce] = generateBC_CZM_BilinearV2(nodeBou, node(:, 2:3), uP);
BC = setBC(fixNode, nodeForce, size(node, 1) * 2);

% uv should satisfy the boundary condition
uv(BC.DirchletDOF) = BC.Dirichlet;
end

function [KCZM, Fcv] = CohesiveMatrix(uv, sigma_c, lamda_cr, delta_c, elem, node, Para)
% Preparation for cohesive element
elemCoh = elem(elem(:, 2) == 2, 3 : end);

uv = reshape(uv, 2, [])';
GaussInfo = shapeFunc_valueDeriv_CZM(elemCoh, node(:, 2:3), uv, Para);
[KCZM, Fcv] = globalK2D_CZM_Bilinear(Para, elemCoh, GaussInfo, uv, 'bilinear', sigma_c, lamda_cr, delta_c);
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