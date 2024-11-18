function [GaussInfo] = initial_DMax(elem, Para)
% -------------------------------------------------------------------
% Precomputing the max seperation displacements at Gauss Pts
% For cohesive zone model element.
% ---------------------------------------------------------------------

GaussInfo = struct;
numEle = size(elem, 1); % num of ele
dim = Para.ndim;

NGPs = [2, 1]; % full integration
[gp, wgt] = gauss_quadrature(NGPs(1), NGPs(2));

ValMax = cell(numEle, 1); 

for ei = 1 : numEle
    ValMax{ei} = zeros(size(gp,1) * dim, 1);
end

GaussInfo.ValMax = ValMax;


end
