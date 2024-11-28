function DMax = renew_DMax(Para, elem, elemIdxCoh, GaussInfo, u, DMax)
numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele
ndim = 2;
numEDofs = numEleNd * ndim;
NNd = Para.NNd;

KVals = zeros(numEDofs^2, numEle); % store the stiff matrix
FVals = zeros(numEDofs, numEle); % store the node force vector

for ei = 1 : numEle
    if ismember(ei, elemIdxCoh)
        Ucoord = u(elem(ei,:),:);
        Ucoord = reshape(Ucoord', [], 1);

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
            Delta_n = DELTA(2);
            Delta_t = DELTA(1);

            if abs(Delta_n) > DMax.ValMax{ei}(gpti * 2 - 1)
                DMax.ValMax{ei}(gpti * 2 - 1) = abs(Delta_n);
            end
            if abs(Delta_t) > DMax.ValMax{ei}(gpti * 2)
                DMax.ValMax{ei}(gpti * 2) = abs(Delta_t);
            end
        end

    end


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