function [K, F] = globalK2D_CZM_Exponential(Para, elem, GaussInfo, u, czmType, varargin)
numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele
ndim = 2;
numEDofs = numEleNd * ndim;
NNd = Para.NNd;
if strcmp(czmType, 'Exp')
    s_c = varargin{1};
    l_cr = varargin{2};
    d_c = varargin{3};
end
KVals = zeros(numEDofs^2, numEle); % store the stiff matrix
FVals = zeros(numEDofs, numEle); % store the node force vector

beta = 1;

for ei = 1 : numEle
    elei = elem(ei,:);
    Ucoord = u(elem(ei,:),:);
    Ucoord = reshape(Ucoord', [], 1);

    Ke = zeros(numEDofs); % element stiff-u
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
        d_n = DELTA(2) + 1e-16;
        d_s = DELTA(1) + 1e-16;
        delta = sqrt((beta * d_s) ^ 2 + d_n ^ 2);

        % Now we only consider loading
        dPhi1 = exp(1) * s_c * delta / d_c * exp(- delta / d_c);
        % temp1 = dPhi1 / delta
        temp1 =  exp(1) * s_c / d_c * exp(- delta / d_c);
        % dPhi2 = dPhi1 / delta * (1 - delta / d_c);
        dPhi2 = exp(1) * s_c / d_c * exp(- delta / d_c) * (1 - delta / d_c);

        C_ss = temp1 * (beta ^ 2) + (beta ^ 4) * (d_s ^ 2) / (delta ^ 2) * (dPhi2 - temp1);
        C_nn = dPhi1 / delta + (d_n ^ 2) / (delta ^ 2) * (dPhi2 - temp1);
        C_sn = (beta ^ 2) * d_s * d_n / (delta ^ 2) * (dPhi2 - temp1);
        C_ns = C_sn;

        t_s = temp1 * (beta ^ 2) * d_s;
        t_n = temp1 * d_n;

        if d_n < 0
            C_nn = C_nn * 1e8;
            t_n = C_nn * d_n;

% t_n = 1e8 * t_n;
        end

        Tc = [t_s; t_n];
        D = [C_ss, C_sn;
            C_ns, C_nn];

        % compute element stiffness at quadrature point
        Ke = Ke + B' * D * B * JW(gpti);
        Fe = Fe + B' * Tc * JW(gpti);

    end
    KVals(:, ei) = Ke(:);
    FVals(:, ei) = Fe(:);
end

J = repmat(1:numEDofs, numEDofs, 1);
I = J';
El = elem';
Eldofs = [El(:)*2-1 El(:)*2]';
ElConn = reshape(Eldofs(:),numEleNd*2, numEle);
ElConn = ElConn';

ii = ElConn(:, I(:))';
jj = ElConn(:, J(:))';

K = sparse(ii(:), jj(:), KVals(:), NNd * 2, NNd * 2);
K = (K + K')/2;

J = repmat(1:numEDofs, 1, 1);
I = J';
El = elem';
Eldofs = [El(:)*2-1 El(:)*2]';
ElConn = reshape(Eldofs(:), numEleNd*2, numEle);
ElConn = ElConn';

ii = ElConn(:, I(:))';

F = sparse(ii(:), 1, FVals(:), NNd * 2, 1);

end