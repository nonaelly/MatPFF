function [K, indIteration] = globalK2D_CZM_Bilinear(Para, elem, GaussInfo, u, czmType, varargin)
numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele
numEDofs = numEleNd * Para.ndim;
if strcmp(czmType, 'bilinear')
    s_c = varargin{1};
    l_cr = varargin{2};
    d_c = varargin{3};
end
indIteration = 0;
KVals = zeros(numEDofs^2, numEle); % store the stiff matrix

for ei = 1 : numEle
    elei = elem(ei,:);
    Ucoord = u(elem(ei,:),:);
    Ucoord = reshape(Ucoord', [], 1);

    Ke = zeros(numEDofs); % element stiff-u

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
        if l_e < l_cr
            C_ss = s_c/(l_cr*d_c);
            C_nn = C_ss;
            C_sn = 0;
            C_ns = 0;
        elseif l_e >= l_cr && l_e < 1
            indIteration = 1;
            C_ss = - d_c*s_c/(1-l_cr)*((d_s/(l_e*d_c^2))^2)...
                + (1-l_e)*(d_c*s_c)/(1-l_cr)*(1/(l_e*d_c^2)-1/(l_e^3)*(d_s^2/(d_c^4)));
            C_nn = - d_c*s_c/(1-l_cr)*((d_n/(l_e*d_c^2))^2)...
                + (1-l_e)*(d_c*s_c)/(1-l_cr)*(1/(l_e*d_c^2)-1/(l_e^3)*(d_n^2/(d_c^4)));
            C_sn = - d_c*s_c/(1-l_cr)*(d_s/(l_e*d_c^2))*(d_n/(l_e*d_c^2))...
                + (1-l_e)*(d_c*s_c)/(1-l_cr)*(-1/(l_e^3)*(d_s*d_n/(d_c^4)));
            C_ns = C_sn;
        else
            C_ss = 0;
            C_nn = 0;
            C_sn = 0;
            C_ns = 0;
        end
        D = [C_ss, C_sn;
            C_ns, C_nn];
        % compute element stiffness at quadrature point
        Ke = Ke + B' * D * B * JW(gpti);
    end
    KVals(:, ei) = Ke(:);
end

J = repmat(1:numEDofs, numEDofs, 1);
I = J';
El = elem';
Eldofs = [El(:)*2-1 El(:)*2]';
ElConn = reshape(Eldofs(:),numEleNd*2, numEle);
ElConn = ElConn';

ii = ElConn(:, I(:))';
jj = ElConn(:, J(:))';

K = sparse(ii(:), jj(:), KVals(:));
K = (K + K')/2;

end