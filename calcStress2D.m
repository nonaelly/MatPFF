function [Stress] = calcStress2D(GaussInfo, elem, Para, Disp)
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

CV = [ 1  1  1  0  0  0;
       1  1  1  0  0  0;
       1  1  1  0  0  0;
       0  0  0  0  0  0;
       0  0  0  0  0  0;
       0  0  0  0  0  0]; % volumetric stiff tr()I
   
CD = [ 4/3 -2/3 -2/3 0  0  0;
      -2/3  4/3 -2/3 0  0  0;
      -2/3 -2/3  4/3 0  0  0;
       0    0    0   1  0  0;
       0    0    0   0  1  0;
       0    0    0   0  0  1]; % deviatoric stiff 2e
   
E = Para.E; % elastic modulus
nu = Para.nu; % poisson's ratio
K = E/(3*(1-2*nu));
G = E/(2*(1+nu));

D = CV * K + CD * G;

numEleNd  = size(elem, 2);  % 单元结点数
numEle = size(elem, 1); % 单元数
numEDofs = numEleNd * Para.ndim;
StressVals = zeros(numEleNd * 6, numEle);
FatigueVal = zeros(numEleNd, numEle);
for ei = 1 : numEle
    elei = elem(ei,:);
    eleDOFs = reshape([2*elei-1; 2*elei], numEDofs,1);
    EleDisp = Disp(eleDOFs);
    
    % loading FEM information
    dRdxGaussPt = GaussInfo.SpDeriv{ei};
    RGaussPt = GaussInfo.SpVal{ei};
    % JW = GaussInfo.JW{ei};
    GStress = zeros(6, size(dRdxGaussPt,3));
    NMtx = zeros(size(dRdxGaussPt,3), numEleNd);
    
    for gpti = 1 : size(dRdxGaussPt,3)
        
        %Compute derivatives of basis functions w.r.t physical coordinates
        dRdx = dRdxGaussPt( :, :, gpti);
        R = RGaussPt(gpti, :);
        
        B = zeros(6, 2 * numEleNd);
        B(1, 1 : 2 : 2 * numEleNd) = dRdx(1, :);
        B(2, 2 : 2 : 2 * numEleNd) = dRdx(2, :);
        B(4, 1 : 2 : 2 * numEleNd) = dRdx(2, :);
        B(4, 2 : 2 : 2 * numEleNd) = dRdx(1, :);
        
        GStress(:,gpti) = D * B * EleDisp;
        NMtx(gpti,:) = R;
    end
    NodeStress = NMtx \ GStress'; % extrapolation % stress on nodes
    StressVals(:, ei) = NodeStress(:);
end

kk = reshape(elem',[],1);

for ii = 1:6
    vals = StressVals([numEleNd*ii-numEleNd+1 : numEleNd*ii], :);
    S(:,ii) = sparse(kk, 1, vals(:));  % assemble rhs
end
Nrepeats = sparse(kk, 1, ones(numel(vals),1)); % 每个结点的叠加次数
S = full(S);
S = S ./ Nrepeats;

Stress.xx = S(:,1);
Stress.yy = S(:,2);
Stress.zz = S(:,3);
Stress.xy = S(:,4);
Stress.yz = S(:,5);
Stress.xz = S(:,6);

% Calculate vonMises stress
firstTerm = 1 / 2 * ((Stress.xx - Stress.yy) .^ 2 + ...
    (Stress.yy - Stress.zz) .^ 2 + ...
    (Stress.zz - Stress.xx) .^ 2);
secondTerm = 3 * (Stress.xy .^ 2 + Stress.yz .^ 2 + Stress.xz .^ 2);
Stress.vonMises = sqrt(firstTerm + secondTerm);

end