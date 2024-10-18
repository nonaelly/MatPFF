function [GaussInfo] = shapeFunc_valueDeriv_CZM(elem, node, u, Para)
% -------------------------------------------------------------------
% Precomputing the shape function's values and derivatives at Gauss Pts
% For: physical value at GPt (shape function's values * nodal values)
%      strain-disp matrix at GPt (shape function's derivatives)
% For cohesive zone model element.
% ---------------------------------------------------------------------

GaussInfo = struct;
Dof = Para.ndim;

numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele

NGPs = [2, 1]; % full integration
[gp, wgt] = gauss_quadrature(NGPs(1), NGPs(2));
wgt = wgt / 2;

EleShapeDerivPara = cell(numEle, 1);  % Ele: shape function's deriv at GPts
JW = cell(numEle, 1);             % Ele: det(Jacobian) * weight at GPts
Hisy = cell(numEle, 1);           % Ele: reference energy at GPts
H = cell(numEle, 1);           % Ele: reference energy at GPts
% JA = cell(numEle, 1);             % Ele: det(Jacobian) at GPts

for ei = 1 : numEle
    Ucoord = u(elem(ei,:),:);
    Ndcoord = node(elem(ei,:),:);

    JacWeight = zeros(size(gp,1),1);
    HMatrix = zeros(size(gp,1),2);
    HisyGaussPt = zeros(size(gp,1),1); % predefine crack with nodal value, HisyGaussPt=0

    delX = (Ucoord(2, 1) + Ucoord(3, 1) + Ndcoord(2, 1) + Ndcoord(3, 1))/2 - ...
        (Ucoord(1, 1) + Ucoord(4, 1) + Ndcoord(1, 1) + Ndcoord(4, 1))/2;
    delY = (Ucoord(2, 2) + Ucoord(3, 2) + Ndcoord(2, 2) + Ndcoord(3, 2))/2 - ...
        (Ucoord(1, 2) + Ucoord(4, 2) + Ndcoord(1, 2) + Ndcoord(4, 2))/2;

    J11 = delX / 2;
    J12 = delY / 2;
    dxdxi = [J11, J12]; % Jacobi matrix
    J = sqrt(J11^2 + J12^2);  % jacobian
    dxdxiGaussPt = dxdxi; %

    for gpti = 1 : size(gp,1)
        Htemp = [(1-gp(gpti))/2, (1+gp(gpti))/2];
        W = wgt(gpti);        % weigths for each integration point
        JacWeight(gpti) = J * W;
        HMatrix(gpti, :) = Htemp;
        %         Jacob(gpti) = J;
    end

    EleShapeDerivPara{ei} = dxdxiGaussPt;
    JW{ei} = JacWeight;
    Hisy{ei} = HisyGaussPt; % all zeros
    H{ei} = HMatrix; % all zeros
    %     JA{ei} = Jacob;
end

GaussInfo.SpDerivPara = EleShapeDerivPara;
GaussInfo.JW = JW;
GaussInfo.Hisy = Hisy;
GaussInfo.H = H;

end
