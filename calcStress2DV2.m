function [Stress] = calcStress2DV2(GaussInfo, elem, Para, Disp)
% -------------------------------------------------------------------
% Calculate two dimensional stiffness matrices
% ------------------------------------------------------------------
% Input:
%       Mesh: Mesh structure
%       NURBS: NURBS structure
%       E: Young's modulus
%       nu: Poisson's ratio
% -------------------------------------------------------------------
% Output:
%       KVals: matrix stores local stiffness matrices
%           (size(KVals) = [(NEN * 2) ^ 2, NEL])
% -------------------------------------------------------------------
% Consider B-Bar

E = Para.E; % elastic modulus
nu = Para.nu; % poisson's ratio
K = E/(3*(1-2*nu));
G = E/(2*(1+nu));

D = E*(1-nu)/((1+nu)*(1-2*nu))*[1,nu/(1-nu),0;nu/(1-nu),1,0;0,0,(1-2*nu)/(2*(1-nu))];

numEleNd  = size(elem, 2);  % 单元结点数
numEle = size(elem, 1); % 单元数
numEDofs = numEleNd * Para.ndim;
StressVals = zeros(numEleNd * 3, numEle);
for ei = 1 : numEle
    elei = elem(ei,:);
    eleDOFs = reshape([2*elei-1; 2*elei], numEDofs,1);
    EleDisp = Disp(eleDOFs);

    % loading FEM information
    dRdxGaussPt = GaussInfo.SpDeriv{ei};
    RGaussPt = GaussInfo.SpVal{ei};
    % JW = GaussInfo.JW{ei};
    GStress = zeros(3, size(dRdxGaussPt,3));
    NMtx = zeros(size(dRdxGaussPt,3), numEleNd);

    dRdx_0 = GaussInfo.EleShapeDerivBar{ei};
    BBarDil = zeros(3, 2 * numEleNd);
    BBarDil(1:2, 1 : 2 : 2 * numEleNd) = 1/2 * [dRdx_0(1, :); dRdx_0(1, :)];
    BBarDil(1:2, 2 : 2 : 2 * numEleNd) = 1/2 * [dRdx_0(2, :); dRdx_0(2, :)];

    for gpti = 1 : size(dRdxGaussPt,3)

        %Compute derivatives of basis functions w.r.t physical coordinates
        dRdx = dRdxGaussPt( :, :, gpti);
        R = RGaussPt(gpti, :);

        B = zeros(3, 2 * numEleNd);
        B(1, 1 : 2 : 2 * numEleNd) = dRdx(1, :);
        B(2, 2 : 2 : 2 * numEleNd) = dRdx(2, :);
        B(3, 1 : 2 : 2 * numEleNd) = dRdx(2, :);
        B(3, 2 : 2 : 2 * numEleNd) = dRdx(1, :);

        BDil = zeros(3, 2 * numEleNd);
        BDil(1:2, 1 : 2 : 2 * numEleNd) = 1/2 * repmat(dRdx(1, :), 2, 1);
        BDil(1:2, 2 : 2 : 2 * numEleNd) = 1/2 * repmat(dRdx(2, :), 2, 1);

        BBar = B - BDil + BBarDil;

        %         GStress(:,gpti) = D * B * EleDisp;
        GStress(:,gpti) = D * BBar * EleDisp;
        NMtx(gpti,:) = R;
    end
    NodeStress = NMtx \ GStress'; % extrapolation % stress on nodes
    StressVals(:, ei) = NodeStress(:);
end

kk = reshape(elem',[],1);

for ii = 1:3
    vals = StressVals([numEleNd*ii-numEleNd+1 : numEleNd*ii], :);
    S(:,ii) = sparse(kk, 1, vals(:));  % assemble rhs
end
Nrepeats = sparse(kk, 1, ones(numel(vals),1)); % 每个结点的叠加次数
S = full(S);
S = S ./ Nrepeats;

Stress.xx = S(:,1);
Stress.yy = S(:,2);
Stress.zz = nu * (S(:,1) + S(:,2));
Stress.xy = S(:,3);
Stress.yz = zeros(size(S, 1), 1);
Stress.xz = zeros(size(S, 1), 1);

% Calculate vonMises stress
firstTerm = 1 / 2 * ((Stress.xx - Stress.yy) .^ 2 + ...
    (Stress.yy - Stress.zz) .^ 2 + ...
    (Stress.zz - Stress.xx) .^ 2);
secondTerm = 3 * (Stress.xy .^ 2 + Stress.yz .^ 2 + Stress.xz .^ 2);
Stress.vonMises = sqrt(firstTerm + secondTerm);

end