% Elastic problem.
% One cohesive element

clear; close all; clc
addpath("Func\")
%%  ***  Reas Ansys Mesh  ***
a = 1;
b = 0;
O = [0,0];
numX = 1;
numY = 1;
% numR = 1;
% numTheta = 1;
[node, elem, nodeBou, elemBou] = generateMeshFEM('rect', a, b, O, numX, numY);
sumNode = size(node,1);

%% ***  Material para  *** (Ambati's Paper)
% Para.PFModel = 2; % 1-AT2; 2-AT1
Para.ndim = 2; % dim
Para.isStress = 2;  % 1 - plane stress, 2 - plane strain
Para.E = 25714; % Young's Modulus based on (N/mm2)
Para.nu = 0.2857; % Poisson's Ratio
Para.lambda = Para.E*Para.nu/((1+Para.nu)*(1-2*Para.nu)); % Lame Constant
Para.mu = Para.E/(2*(1+Para.nu)); % Lame Constant
% Para.Gc = 0.089; % Critical energy release for unstable crack (Gc, N/mm)
% Para.Len = 3; %

Para.NNd = size(node,1); % number of nodes

elem(:,1:2) = [];
node(:,1)   = [];
node = node(:, 1 : Para.ndim);
u = zeros(size(node));
[GaussInfo] = shapeFunc_valueDeriv_CZM(elem, node, u, Para);
%% Cohesive zone model

sMax = 1;
delta = 0.1;
alphaP = 2;
KCZM = globalK2D_CZM(Para, elem, GaussInfo, u, 1, sMax, delta, alphaP);

%% Elastic problem
K = globalK2D(Para, elem, GaussInfo);

% Set boundary condition.
boundary = {'p', 'dx', 'dy', 'f'};
pressure = 1;
[fixNode, nodeForce] = generateBC(boundary, nodeBou, elemBou, elem, node, pressure);

BC = setBC(fixNode, nodeForce, sumNode*2);

% solve
Disp = solveBC(BC, K);
Stress = calcStress2DV2(GaussInfo, elem, Para, Disp);

Disp = reshape(Disp,Para.ndim,[]);
Disp = Disp';

r = b;
uOut = (1+Para.nu)*pressure*a^2/Para.E/(b^2-a^2)*((1-2*Para.nu)*r + b^2/r);

%% Plot
% ux
figure
axis equal;
PlotContour(node,elem,Disp(:,1),'ux',1);
axis off;
% uy
figure
axis equal;
PlotContour(node,elem,Disp(:,2),'uy',1);
axis off;
% mises
figure
axis equal;
PlotContour(node,elem,Stress.vonMises,'mises',1);
axis off;

%% Viscoelastic problem
fprintf(1,'sort the global K1, K2 for vis\n')
K1 = globalK2D(Para, elem, GaussInfo);
[K2, K3] = globalKVis2D(Para, elem, GaussInfo);

% Solve viscoelastic process

t0=0;   % start time
tEnd=300;  % total time
dt=0.1; % time step

E = Para.E; % elastic modulus
nu = Para.nu; % poisson's ratio
K0 = E/(3*(1-2*nu));
G0 = E/(2*(1+nu));
tauG = 1/2.3979;  % relaxation times
tauK = 1/2.3979;  % relaxation times
visG = 0.99; % nomalized shear modulus
visK = 0;    % nomalized volume modulus

[qG_jk, qG_jk1] = deal(zeros(numel(node), length(visG)));
[qK_jk, qK_jk1] = deal(zeros(numel(node), length(visK)));
[uk1, uk2, Fvis] = deal(zeros(numel(node), 1));
uTest = zeros(numel(node), (tEnd-t0)/dt+1);
t = t0;
k = 1;
while k <= (tEnd-t0)/dt+1
    g_dt = (1 - visG.*(1-exp(-dt./tauG))) ;
    k_dt = (1 - visK.*(1-exp(-dt./tauK)));

    if k == 1
        KA = K1;
    else
        KA = K1 + G0*(1 - g_dt)*K2 + K0*(1 - k_dt)*K3;
        % u (k-1) = u (k)
        uk2 = uk1;
        uk1 = u;
        % u* (k-2) = 1/2*(u (k-2) +  u (k-1))
        uStar_k2 = (uk2 + uk1)/2;

        % q_j,K(k) = exp(-dt/xi)*((1-exp(-dt/xi) * u* + q_j,K(k-1)))
        qG_jk1 = qG_jk;
        qK_jk1 = qK_jk;
        if k > 2
            qG_jk = exp(-dt./tauG).*((1-exp(-dt./tauG)).*uStar_k2+qG_jk1);
            qK_jk = exp(-dt./tauK).*((1-exp(-dt./tauK)).*uStar_k2+qK_jk1);
        end
    end

    feM1 = - 2*K2*(G0*sum(visG.*qG_jk, 2) + 1/2*G0*(1 - g_dt)*uk1);
    feM2 = - K3*(K0*sum(visK.*qK_jk, 2) + 1/2*K0*(1 - k_dt)*uk1);
    %     force = feT + feM + feM1 + feM2 + feM3;
    Fvis = feM1 + feM2;

    u = solveBC(BC, KA, Fvis);
    uTest(:,k) = u;

    t = t + dt;
    k = k + 1;
end

uAna = CylinderAnalysis([0, b], a, b, -pressure, (1-visG)*G0, visG*G0, Para.nu, t0:dt:tEnd, tauG);

figure
n = 100;
plot(t0:n*dt:tEnd, uAna(1:n:end, 2), 'bo')
hold on
plot(t0:dt:tEnd, uTest(2,:), 'r-')
