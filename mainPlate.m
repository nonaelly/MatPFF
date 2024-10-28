% Elastic problem.
% two different elastic plate. 
% Compared with an analytical solution.

clear; close all; clc
addpath("Func\")
%%  ***  Reas Ansys Mesh  ***
L = 1;
a = L*[1, 1, 1];
b = [0.5, 0.5, 0];
O = [0,0; 0,0.5; 0,0.5];
numX = 4*ones(3, 1);
numY = [15, 15, 1];
[node, elem, nodeBou, elemBou] = deal(cell(3, 1));
sumNode = zeros(3, 1);
for i = 1:3
    [node{i}, elem{i}, nodeBou{i}, elemBou{i}] = generateMeshFEM('rect', a(i), b(i), O(i,:), numX(i), numY(i));
    sumNode(i) = size(node{i},1);
end

%% ***  Material para  *** (Ambati's Paper)
E = [100, 200, 0];
[Para, GaussInfo] = deal(cell(3, 1));
for i = 1:3
    Para{i}.ndim = 2; % dim
    Para{i}.isStress = 2;  % 1 - plane stress, 2 - plane strain
    Para{i}.E = E(i); % Young's Modulus based on (N/mm2)
    Para{i}.nu = 0; % Poisson's Ratio
    Para{i}.lambda = Para{i}.E*Para{i}.nu/((1+Para{i}.nu)*(1-2*Para{i}.nu)); % Lame Constant
    Para{i}.mu = Para{i}.E/(2*(1+Para{i}.nu)); % Lame Constant
    Para{i}.NNd = size(node{i},1); % number of nodes

    elem{i}(:,1:2) = [];
    node{i}(:,1)   = [];
    node{i} = node{i}(:, 1 : Para{i}.ndim);
    if i == 3
        u = zeros(size(node{i}));
        [GaussInfo{i}] = shapeFunc_valueDeriv_CZM(elem{i}, node{i}, u, Para{i});
    else
        [GaussInfo{i}] = shapeFunc_valueDeriv(elem{i}, node{i}, Para{i});
    end

end

%% Elastic problem with Cohesive zone model
gc_I = 0.1;
tu = 1;
delatn = gc_I / (tu * exp(1));
KCZM = globalK2D_CZM(Para{3}, elem{3}, GaussInfo{3}, u, 'exa', gc_I, delatn);

K = cell(2, 1);
for i = 1 : 2
    K{i} = globalK2D(Para{i}, elem{i}, GaussInfo{i});
end

% Assemble the total K.
Ktotal = blkdiag(K{1}, K{2});
nodeTotal = [node{1}; node{2}];
idx =  size(K{1}, 1) - size(KCZM, 1)/2 + 1 : size(K{1}, 1) + size(KCZM, 1)/2;
Ktotal(idx, idx) = Ktotal(idx, idx) + KCZM;

Fc_0 = nodeForce_CZM(Para, elem, GaussInfo, u_c0, 'bilinear', sigma_c, lamda_cr, delta_c);

% Set boundary condition.
ubar = 1e-4;
[fixNode, nodeForce] = generateBC_CZM(nodeBou, node, ubar);

BC = setBC(fixNode, nodeForce, sum(sumNode(1:2))*2);

% solve
Disp = solveBC(BC, Ktotal);
% Stress = calcStress2DV2(GaussInfo, elem, Para, Disp);

Disp = reshape(Disp,Para{1}.ndim,[]);
Disp = Disp';

figure
% x = 0
idx_x0 = 1:numX(1)+1:size(nodeTotal, 1);
vn = delatn;
n = 20;
y = zeros(n+1, 1);
for i = 1 : n+1
    y(i) = (i-1) / n * L;
    if y(i) >= 1/2*L
        DispAna(i) = ubar - ubar / ((L/2)*(1+E(2)/E(1)) + E(2)*(vn^2)/gc_I) * (L-y(i));
    else
        DispAna(i) = E(2) * ubar / ((L/2)*(E(1)+E(2)) + E(1)*E(2)*(vn^2)/gc_I) * y(i);
    end
end
plot(nodeTotal(idx_x0, 2), Disp(idx_x0, 2));
hold on
plot(y, DispAna, 'o')

figure
plot(nodeTotal(idx_x0, 2), Disp(idx_x0, 1));

%% Plot
% uy
figure
axis equal;
PlotContour(node{1},elem{1},Disp(1:size(node{1}, 1), 2),'uy',1);
hold on
PlotContour(node{2},elem{2},Disp(end - size(node{2}, 1) + 1 : end, 2),'uy',1);
axis off;

