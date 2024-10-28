function F = nodeForce_CZM(Para, elem, GaussInfo, u, czmType, varargin)
numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele
numEDofs = numEleNd * Para.ndim;
if strcmp(czmType, 'bilinear')
    s_c = varargin{1};
    l_cr = varargin{2};
    d_c = varargin{3};
end
FVals = zeros(numEDofs, numEle); % store the node force vector

for ei = 1 : numEle
    elei = elem(ei,:);
    Ucoord = u(elem(ei,:),:);
    Ucoord = reshape(Ucoord', [], 1);

    Fe = zeros(numEDofs, 1); % element node force

    % loading FEM information
    dxdxi = GaussInfo.SpDerivPara{ei};
    JW = GaussInfo.JW{ei};
    H1H2 = GaussInfo.H{ei};

    J11 = dxdxi(1);
    J12 = dxdxi(2);
    delX = J11 * 2;
    delY = J12 * 2;
    cosTheta = delX / sqrt(delX^2 + delY^2);
    sinTheta = delY / sqrt(delX^2 + delY^2);
    LAMADA = [cosTheta sinTheta; -sinTheta cosTheta];
    R = [];
    for i=1:4
        R = blkdiag(R, LAMADA);
    end

    % Δ1, ... ,Δ4
    % For 0 thick element.
    L = zeros(4, 8);
    temp = [7, 8, 5, 6];
    for i = 1:4
        L(i, i) = -1;
        L(i, temp(i)) = 1;
    end

    for gpti = 1 : size(JW,1)

        H1 = H1H2(gpti, 1);
        H2 = H1H2(gpti, 2);
        H = [H1, 0, H2, 0; 0, H1, 0, H2];

        B = H * L * R;

        DELTA = B * Ucoord;
        d_n = DELTA(2);
        d_s = DELTA(1);
        l_e = sqrt((d_n/d_c)^2 + (d_s/d_c)^2);

        if l_e <= l_cr
            t_s = s_c / l_cr * (d_s / d_c);
            t_n = s_c / l_cr * (d_n / d_c);
            Tc = [t_s; t_n];
        elseif l_e > l_cr && l_e <= 1
            t_s = s_c * (1 - l_e) / (1 - l_cr) / l_e * (d_s / d_c);
            t_n = s_c * (1 - l_e) / (1 - l_cr) / l_e * (d_n / d_c);
            Tc = [t_s; t_n];
        else
            Tc = zeros(2, 1);
        end

        % compute element stiffness at quadrature point
        Fe = Fe + B' * Tc * JW(gpti);
    end
    FVals(:, ei) = Fe(:);
end

J = repmat(1:numEDofs, 1, 1);
I = J';
El = elem';
Eldofs = [El(:)*2-1 El(:)*2]';
ElConn = reshape(Eldofs(:), numEleNd*2, numEle);
ElConn = ElConn';

ii = ElConn(:, I(:))';

F = sparse(ii(:), 1, FVals(:));

end