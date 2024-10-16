function Disp = solveBC(BC, K, Fvis)
if nargin < 3
    F = BC.RHS;
else
    F = BC.RHS + Fvis;
end
% Big number
max_abs_value = max(max(abs(K)));
n = floor(log10(max_abs_value));
bigN = 10^(9+n);
for i = 1 : size(BC.DirchletDOF, 1)
    ind = BC.DirchletDOF(i);
    F(ind) = BC.Dirichlet(i) * K(ind, ind) * bigN;
    K(ind, ind) = K(ind, ind) * bigN;
end
Disp = K\F;
end