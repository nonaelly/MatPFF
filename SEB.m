% single-edge notched beam (SE(B)) test
% <A bilinear cohesive zone model tailored for fracture of asphalt concrete
%   considering viscoelastic bulk material>

clear; close all; clc
addpath("Func\")
%%  ***  Reas Ansys Mesh  ***
% dx = 2;
% numY = [19, 81/1.5];
dx = 4;
numY = [2, 3];
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

numStep = 10;
un = zeros(2*Para.NNd, numStep + 1);
for n = 1 : numStep
    % a starting value for the unknown must be chosen;
    % usually the solution dn from the last time step (n) is selected
    % step of the iteration: v = 0
    uv = un(n, :);
    % Residual
    Rv = K * uv + Fcv - BC.RHS + bigN * uBC;
    % Tanget matrix
    dRduv = K + KCZM + bigN * KBC;
    idxIter = 0;

    while max(abs(Rv)) > tole && idxIter < 100
        % step of the iteration: v = 1, 2, ...
        duv = -dRduv\Rv;
        
        % step of the iteration: v + 1, renew u, R, dRdu
        uv = uv + duv;
        % Residual
        Rv = K * uv + Fcv - BC.RHS + bigN * uBC;
        % Tanget matrix
        dRduv = K + KCZM + bigN * KBC;

        idxIter = idxIter + 1;
    end

    un(n + 1, :) = uv;
    fprintf('Step %d', num2str(n), 'iteration number %d', num2str(idxIter))

end
