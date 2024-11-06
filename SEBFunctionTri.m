function [CMOD, P, u, nodeEla, elemEla] = SEBFunctionTri(YourModel, CZMType, sigma_c, G_c, duP, maxCMOD, varargin)
% single-edge notched beam (SE(B)) test

%%  ***  Read Abaqus Mesh  ***
parts = loadinp(YourModel);

% Copy the half model
rightNode = parts.Nodes;
leftNode = [- rightNode(:, 1), rightNode(:, 2)];

node = [rightNode; leftNode];
node = [(1 : size(node, 1))', node];

% Find cohesive node
coNodes = node((rightNode(:, 1) == 0) & (rightNode(:, 2) >= 19) & (rightNode(:, 2) <= 100), :);
coNodes = sortrows(coNodes, 3);

% Cohesive element
coElme = [];
for i = 1 : size(coNodes, 1) - 1
    [coElme] = [coElme; coNodes(i, 1), coNodes(i + 1, 1), ...
        coNodes(i + 1, 1) + size(rightNode, 1), coNodes(i, 1) + size(rightNode, 1)];
end

rightElem = [];
for i = 1 : size(parts.Elements, 1)
    [rightElem] = [rightElem; parts.Elements(i).Connectivity, parts.Elements(i).Connectivity(end)];
end
leftElem = rightElem + size(rightNode, 1);
leftElem(:, [2, 4]) = leftElem(:, [4, 2]);
[elem] = [rightElem; leftElem];
[elem] = [(1 : size(elem, 1))', ones(size(elem, 1), 1), elem];
[elem] = [elem; (1 + size(elem, 1) : size(elem, 1) + size(coElme, 1))', 2 * ones(size(coElme, 1), 1), coElme];

elemEla = elem(elem(:, 2) == 1, 3 : end);
elemCoh = elem(elem(:, 2) == 2, 3 : end);
nodeEla = node(1 : max(max(elemEla)), 2:3);

% % Plot
% figure(1);
% hold on;
% for i = 1:size(elem, 1)
%     nodes = node(elem(i, 3 : end), 2:3);
%     fill(nodes(:, 1), nodes(:, 2), 'w', 'EdgeColor', 'k');
%     if elem(i, 2) == 2 % CZM
%         plot(node(elem(i, 3:end), 2), node(elem(i, 3:end), 3), 'ro', 'MarkerFaceColor', 'r');
%     end
% end
% 
% title('single-edge notched beam (SE(B)) test');
% xlabel('X');
% ylabel('Y');
% axis equal;

node = node(:, 2 : 3);
elem = elem(:, 3 : end);

%% ***  Material para  *** (Ambati's Paper)
E = 14.2e3;
Para.ndim = 2; % dim
% Para.isStress = 2;  % 1 - plane stress, 2 - plane strain
Para.isStress = 1;  % 1 - plane stress, 2 - plane strain
Para.E = E; % Young's Modulus based on (N/mm2)
Para.nu = 0.35; % Poisson's Ratio
Para.lambda = Para.E*Para.nu/((1+Para.nu)*(1-2*Para.nu)); % Lame Constant
Para.mu = Para.E/(2*(1+Para.nu)); % Lame Constant
Para.NNd = size(node,1); % number of nodes

%% Elastic problem with Cohesive zone model
switch CZMType
    case 'bilinear'
        delta_c = 2*G_c/sigma_c*1e-3; % mm
    case 'Exp'
        delta_c = G_c/sigma_c/exp(1)*1e-3; % mm
end
lamda_cr = varargin{1};

% Generate elastic stiffness matrix
GaussInfo = deal(cell(2, 1));
[GaussInfo{1}] = shapeFunc_valueDeriv_tri(elemEla, node, Para);
K = globalK2D(Para, elemEla, GaussInfo{1});

% Newton-Raphson method
tole = 1e-6;

uP = 0;
un = zeros(2*Para.NNd, 1);
CMOD = 0;
P = 0;
idxCMOD = find(ismember(node, [0, 0], 'rows'));
n = 1;
u = [];
while abs(CMOD(n)) <= maxCMOD
    % Boundary condition
    uP = uP - duP;

    % A starting value for the unknown must be chosen;
    % usually the solution u_n from the last time step (n) is selected
    % modify u_n to satisfy the boundary conditions.

    % Step of the iteration: v = 0
    [BC, uv] = startingValue(un, uP, node);

    % Calculate the cohesive nodal forces and the stiffness matrix
    [KCZM, Fcv] = CohesiveMatrix(uv, sigma_c, delta_c, elemCoh, node, Para, CZMType, lamda_cr);

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
        [KCZM, Fcv] = CohesiveMatrix(uv, sigma_c, delta_c, elemCoh, node, Para, CZMType, lamda_cr);

        % Residual
        Rv = K * uv + Fcv - BC.RHS;
        % Tanget matrix
        dRduv = K + KCZM;

        tempP = sum(Rv(BC.DirchletDOF(abs(BC.Dirichlet) > 0))) * 75 / 1e3;
        %         tempP = (Rv(BC.DirchletDOF(1)));

        % Modify Rv and dRduv to satisfy the boundary conditions.
        [Rv, dRduv] = ApplyBC(BC, Rv, dRduv);

        idxIter = idxIter + 1;
    end

    un = uv;
    [u] = [u, un];
    n = n + 1;

    CMOD(n) = uv(idxCMOD(1) * 2 - 1) - uv(idxCMOD(2) * 2 - 1);
    P(n) = tempP;

    fprintf('Step %d, iteration number %d, uP = %f, P = %f, CMOD= %f\n', n, idxIter, uP, tempP, CMOD(n));
end

end

%% Sub function
function [BC, uv] = startingValue(un, uP, node)
uv = un;

% Boundary condition
fixNode = [];
nodeForce = [];

% (162, 0) and (- 162, 0)
m = find(ismember(node, [162, 0], 'rows'));
[fixNode] = [fixNode; m, 1, 0; m, 2, 0];
m = find(ismember(node, [-162, 0], 'rows'));
[fixNode] = [fixNode; m, 2, 0];

% (0, 100) P
m = find(ismember(node, [0, 100], 'rows'));
[fixNode] = [fixNode; m(1), 2, uP; m(2), 2, uP];

fixNode = sortrows(fixNode, 1);

BC = setBC(fixNode, nodeForce, size(node, 1) * 2);

% uv should satisfy the boundary condition
uv(BC.DirchletDOF) = BC.Dirichlet;
end

function [KCZM, Fcv] = CohesiveMatrix(uv, sigma_c, delta_c, elemCoh, node, Para, CZMType, lamda_cr)

uv = reshape(uv, 2, [])';
GaussInfo = shapeFunc_valueDeriv_CZM(elemCoh, node, uv, Para);
switch CZMType
    case 'bilinear'
        [KCZM, Fcv] = globalK2D_CZM_Bilinear(Para, elemCoh, GaussInfo, uv, CZMType, sigma_c, lamda_cr, delta_c);
    case 'Exp'
        [KCZM, Fcv] = globalK2D_CZM_Exponential(Para, elemCoh, GaussInfo, uv, CZMType, sigma_c, lamda_cr, delta_c);
end

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