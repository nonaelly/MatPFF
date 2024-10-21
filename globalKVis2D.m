function [K2, K3] = globalKVis2D(Para, elem, GaussInfo)
isStress = Para.isStress;  % 1 - plane stress, 2 - plane strain
E = Para.E; % Young's Modulus based on (N/mm2)
nu = Para.nu; % Poisson's Ratio

D2=[-4/3  2/3 0;
    2/3 -4/3 0;
    0    0  -1];
D3=[-1 -1 0;
    -1 -1 0;
    0  0 0];

numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele
numEDofs = numEleNd * Para.ndim;

K2Vals = zeros(numEDofs^2, numEle); % store the stiff matrix
K3Vals = zeros(numEDofs^2, numEle); % store the stiff matrix

for ei = 1 : numEle
    elei = elem(ei,:);
    K2e = zeros(numEDofs); % element stiff-u
    K3e = zeros(numEDofs); % element stiff-u

    % loading FEM information
    dRdxGaussPt = GaussInfo.SpDeriv{ei};
    RGaussPt = GaussInfo.SpVal{ei};
    JW = GaussInfo.JW{ei};

    dRdx_0 = GaussInfo.EleShapeDerivBar{ei};
    R0 = GaussInfo.EleShapeValBar{ei};
    BBarDil = zeros(3, 2 * numEleNd);
    BBarDil(1:2, 1 : 2 : 2 * numEleNd) = 1/2 * [dRdx_0(1, :); dRdx_0(1, :)];
    BBarDil(1:2, 2 : 2 : 2 * numEleNd) = 1/2 * [dRdx_0(2, :); dRdx_0(2, :)];

    %     % BBar = ∫BdΩ/∫dΩ
    %     inteB = zeros(2, numEleNd);
    %     S = 0;
    %     for gptj = 1 : size(dRdxGaussPt,3)
    %         dRdx = dRdxGaussPt( :, :, gptj);
    %         inteB = inteB + dRdx * JW(gptj) * 1;
    %         S = S + JW(gptj) * 1;
    %     end
    %     BBarDil2 = zeros(3, 2 * numEleNd);
    %     BBarDil2(1:2, 1 : 2 : 2 * numEleNd) = 1/2 /S * repmat(inteB(1, :), 2, 1);
    %     BBarDil2(1:2, 2 : 2 : 2 * numEleNd) = 1/2 /S * repmat(inteB(2, :), 2, 1);

    for gpti = 1 : size(dRdxGaussPt,3)

        %Compute derivatives of basis functions w.r.t physical coordinates
        dRdx = dRdxGaussPt( :, :, gpti);
        R = RGaussPt(gpti, :);

        % %      _                                             _
        % %     | N_{1, x} 0         ... N_{m, x} 0         ...|
        % %  B =| 0        N_{1, y}  ... 0        N_{m, y}  ...|
        % %     | N_{1, y} N_{1, x}  ... N_{m, y} N_{m, x}  ...|
        % %      -                                             -
        % % \sigma = [\sigma_{xx} \sigma_{yy} \sigma_{xy}]
        B = zeros(3, 2 * numEleNd);
        BDil = zeros(3, 2 * numEleNd);
        B(1, 1 : 2 : 2 * numEleNd) = dRdx(1, :);
        B(2, 2 : 2 : 2 * numEleNd) = dRdx(2, :);

        B(3, 1 : 2 : 2 * numEleNd) = dRdx(2, :);
        B(3, 2 : 2 : 2 * numEleNd) = dRdx(1, :);

        BDil(1:2, 1 : 2 : 2 * numEleNd) = 1/2 * repmat(dRdx(1, :), 2, 1);
        BDil(1:2, 2 : 2 : 2 * numEleNd) = 1/2 * repmat(dRdx(2, :), 2, 1);

        BBar = B - BDil + BBarDil;
%         BBar = B;
        %         BBar = B - BDil + BBarDil2;

        % compute element stiffness at quadrature point
        % h = 1
        %                 Ke = Ke + B' * D * B * JW(gpti) * 1;
        K2e = K2e + 1 / 2 * BBar' * D2 * BBar * JW(gpti) * 1;
        K3e = K3e + BBar' * D3 * BBar * JW(gpti) * 1;
    end
    K2Vals(:, ei) = K2e(:);
    K3Vals(:, ei) = K3e(:);
end

J = repmat(1:numEDofs, numEDofs, 1);
I = J';
El = elem';
Eldofs = [El(:)*2-1 El(:)*2]';
ElConn = reshape(Eldofs(:),numEleNd*2, numEle);
ElConn = ElConn';

ii = ElConn(:, I(:))';
jj = ElConn(:, J(:))';

K2 = sparse(ii(:), jj(:), K2Vals(:));
K2 = (K2 + K2')/2;

K3 = sparse(ii(:), jj(:), K3Vals(:));
K3 = (K3 + K3')/2;

end