function [GaussInfo] = initial_StateVaribles(elem, Para)
% -------------------------------------------------------------------
% Precomputing the max seperation displacements at Gauss Pts
% For cohesive zone model element.
% ---------------------------------------------------------------------

GaussInfo = struct;
numEle = size(elem, 1); % num of ele
dim = Para.ndim;

NGPs = [2, 1]; % full integration
[gp, wgt] = gauss_quadrature(NGPs(1), NGPs(2));

valUMax = cell(numEle, 1); 
isFail = cell(numEle, 1); 

for ei = 1 : numEle
    valUMax{ei} = zeros(size(gp,1) * dim, 1);
    isFail{ei} = zeros(size(gp,1), 1);
end

GaussInfo.ValMax = valUMax;
GaussInfo.isFail = isFail;

end
