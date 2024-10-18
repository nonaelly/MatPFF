function K = globalK2D_CZM(Para, elem, GaussInfo, u, czmType, varargin)
isStress = Para.isStress;  % 1 - plane stress, 2 - plane strain
E = Para.E; % Young's Modulus based on (N/mm2)
nu = Para.nu; % Poisson's Ratio


numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele
numEDofs = numEleNd * Para.ndim;
switch czmType
    case 1
        sMax = varargin{1};
        delta = varargin{2};
        alphaP = varargin{3};
        const1 = 27/4*sMax/delta;

end

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
        dn = DELTA(2)/delta;
        dt = DELTA(1)/delta;
        switch czmType
            case 1
                if dn<1
                    D = const1 * [alphaP*(1-2*dn+(dn^2)), alphaP*dt*(-2+2*dn);
                        (1-4*dn+3*(dn^2)) + alphaP*(dt^2), alphaP*2*dt*(dn-1)];
                else
                    D = zeros(2);
                end
        end
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