function [K, F] = globalK2D_CZM_PPR(Para, elem, GaussInfo, u, czmType, varargin)
numEleNd  = size(elem, 2);  % num of ele nodes
numEle = size(elem, 1); % num of ele
ndim = 2;
numEDofs = numEleNd * ndim;
NNd = Para.NNd;
if strcmp(czmType, 'PPR')
    ParaPPR = varargin{1};
    DMax = varargin{2};
end
% PPR parameters
delta_n = ParaPPR(1);
delta_t = ParaPPR(2);
m = ParaPPR(3);
n = ParaPPR(4);
alpha_k = ParaPPR(5);
beta_k = ParaPPR(6);
sigma_max = ParaPPR(7);
tau_max = ParaPPR(8);

% generative parameters
lambda_n = sqrt(m / (alpha_k * (m + alpha_k - 1)));
lambda_t = sqrt(n / (beta_k * (n + beta_k - 1)));
phi_n = sigma_max * delta_t / (alpha_k * lambda_n * (1 - lambda_n)^(alpha_k - 1) * (alpha_k / m + 1) * (alpha_k / m * lambda_n + 1)^(m - 1));
phi_t = tau_max * delta_n / (beta_k * lambda_t * (1 - lambda_t)^(beta_k - 1) * (beta_k / n + 1) * (beta_k / n * lambda_t + 1)^(n - 1));

if phi_n == phi_t
    Gamma_n = - phi_n * (alpha_k / m)^m;
    Gamma_t = (beta_k / n)^n;
else
    Gamma_n = (- phi_n)^(max(phi_n - phi_t, 0) / (phi_n - phi_t)) * (alpha_k / m)^m;
    Gamma_t = (- phi_t)^(max(phi_t - phi_n, 0) / (phi_t - phi_n)) * (beta_k / n)^n;
end
alpha_c = - Gamma_n / (delta_n^2) * ((m / alpha_k) ^ (m - 1)) * (alpha_k + m) ...
    * (Gamma_t * (n / beta_k) ^ n + max(phi_t - phi_n, 0));
% unload / reload parameter
alpha_r = 1;
beta_r = 1;

KVals = zeros(numEDofs^2, numEle); % store the stiff matrix
FVals = zeros(numEDofs, numEle); % store the node force vector

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
        Delta_n = DELTA(2) + 1e-16;
        Delta_t = DELTA(1) + 1e-16;

        D_nMax = DMax.ValMax{ei}(gpti * 2 - 1);
        D_tMax = DMax.ValMax{ei}(gpti * 2);
        if Delta_n > delta_n || (Delta_n <= delta_n && abs(Delta_t) > delta_t)
            % interface failure
            D_nn = 1e-15;
            D_nt = 0;
            D_tt = 1e-15;
            D_tn = 0;

            T_n = 1e-15;
            T_t = 1e-15;
        else
            % normal state judge
            if Delta_n >= D_nMax
                % soften
                D_nn = Gamma_n / delta_n^2 * ((m^2 - m) * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^(m - 2) ...
                    + (alpha_k^2 - alpha_k) * (1 - Delta_n / delta_n)^(alpha_k - 2) * (m / alpha_k + Delta_n / delta_n)^m ...
                    - 2 * alpha_k * m * (1 - Delta_n / delta_n)^(alpha_k - 1) * (m / alpha_k + Delta_n / delta_n)^(m - 1))...
                    * (Gamma_t * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^n + max(phi_t - phi_n, 0));

                D_nt = Gamma_n * Gamma_t / (delta_n * delta_t) * (m * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^(m - 1) ...
                    - alpha_k * (1 - Delta_n / delta_n)^(alpha_k - 1) * (m / alpha_k + Delta_n / delta_n)^m)...
                    * (n * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^(n - 1) ...
                    - beta_k * (1 - abs(Delta_t) / delta_t)^(beta_k - 1) * (n / beta_k + abs(Delta_t) / delta_t)^n) * Delta_t / abs(Delta_t);

                T_n = Gamma_n / delta_n * (m * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^(m - 1) ...
                    - alpha_k * (1 - Delta_n / delta_n)^(alpha_k - 1) * (m / alpha_k + Delta_n / delta_n)^m)...
                    * (Gamma_t * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^n + max(phi_t - phi_n, 0));

            elseif Delta_n < D_nMax && Delta_n >= 0
                % unload / reload
                D_nn = Gamma_n / delta_n * (m * (1 - D_nMax / delta_n)^alpha_k * (m / alpha_k + D_nMax / delta_n)^(m - 1) ...
                    - alpha_k * (1 - D_nMax / delta_n)^(alpha_k - 1) * (m / alpha_k + D_nMax / delta_n)^m)...
                    * (Gamma_t * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^n + max(phi_t - phi_n, 0))...
                    * alpha_r / D_nMax * (Delta_n / D_nMax)^(alpha_r - 1);

                D_nt = Gamma_n * Gamma_t / (delta_n * delta_t) * (m * (1 - D_nMax / delta_n)^alpha_k * (m / alpha_k + D_nMax / delta_n)^(m - 1) ...
                    - alpha_k * (1 - D_nMax / delta_n)^(alpha_k - 1) * (m / alpha_k + D_nMax / delta_n)^m)...
                    * (n * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^(n - 1) ...
                    - beta_k * (1 - abs(Delta_t) / delta_t)^(beta_k - 1) * (n / beta_k + abs(Delta_t) / delta_t)^n) * Delta_t / abs(Delta_t) * (Delta_n / D_nMax)^alpha_r;

                T_n = Gamma_n / delta_n * (m * (1 - D_nMax / delta_n)^alpha_k * (m / alpha_k + D_nMax / delta_n)^(m - 1) ...
                    - alpha_k * (1 - D_nMax / delta_n)^(alpha_k - 1) * (m / alpha_k + D_nMax / delta_n)^m)...
                    * (Gamma_t * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^n + max(phi_t - phi_n, 0))...
                    * (Delta_n / D_nMax)^alpha_r;

            elseif Delta_n < 0
                % contact
                D_nn = alpha_c;
                D_nt = 0;
                T_n = D_nn * Delta_n;
            end

            % tangent state judge
            if Delta_n < 0
                % contact (set Delta_n = 0)
                Delta_n = 0;
            end

            if abs(Delta_t) >= D_tMax
                % soften
                D_tt = Gamma_t / delta_t^2 * ((n^2 - n) * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^(n - 2) ...
                    + (beta_k^2 - beta_k) * (1 - abs(Delta_t) / delta_t)^(beta_k - 2) * (n / beta_k + abs(Delta_t) / delta_t)^n ...
                    - 2 * beta_k * n * (1 - abs(Delta_t) / delta_t)^(beta_k - 1) * (n / beta_k + abs(Delta_t) / delta_t)^(n - 1))...
                    * (Gamma_n * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^m + max(phi_n - phi_t, 0));

                D_tn = Gamma_n * Gamma_t / (delta_n * delta_t) * (m * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^(m - 1) ...
                    - alpha_k * (1 - Delta_n / delta_n)^(alpha_k - 1) * (m / alpha_k + Delta_n / delta_n)^m)...
                    * (n * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^(n - 1) ...
                    - beta_k * (1 - abs(Delta_t) / delta_t)^(beta_k - 1) * (n / beta_k + abs(Delta_t) / delta_t)^n) * Delta_t / abs(Delta_t);

                T_t = Gamma_t / delta_t * (n * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^(n - 1) ...
                    - beta_k * (1 - abs(Delta_t) / delta_t)^(beta_k - 1) * (n / beta_k + abs(Delta_t) / delta_t)^n)...
                    * (Gamma_n * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^m + max(phi_n - phi_t, 0)) * Delta_t / abs(Delta_t);

            elseif abs(Delta_t) < D_tMax
                % if abs(Delta_t) < D_tMax
                % unload / reload
                D_tt = Gamma_t / delta_t * (n * (1 - D_tMax / delta_t)^beta_k * (n / beta_k + D_tMax / delta_t)^(n - 1) ...
                    - beta_k * (1 - D_tMax / delta_t)^(beta_k - 1) * (n / beta_k + D_tMax / delta_t)^n)...
                    * (Gamma_n * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^m + max(phi_n - phi_t, 0))...
                    * beta_r / D_tMax * (abs(Delta_t) / D_tMax)^(beta_r - 1);

                D_tn = Gamma_n * Gamma_t / (delta_n * delta_t) * (m * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^(m - 1) ...
                    - alpha_k * (1 - Delta_n / delta_n)^(alpha_k - 1) * (m / alpha_k + Delta_n / delta_n)^m)...
                    * (n * (1 - D_tMax / delta_t)^beta_k * (n / beta_k + D_tMax / delta_t)^(n - 1) ...
                    - beta_k * (1 - D_tMax / delta_t)^(beta_k - 1) * (n / beta_k + D_tMax / delta_t)^n) * Delta_t / abs(Delta_t) * (abs(Delta_t) / D_tMax)^beta_r;

                T_t = Gamma_t / delta_t * (n * (1 - D_tMax / delta_t)^beta_k * (n / beta_k + D_tMax / delta_t)^(n - 1) ...
                    - beta_k * (1 - D_tMax / delta_t)^(beta_k - 1) * (n / beta_k + D_tMax / delta_t)^n)...
                    * (Gamma_n * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^m + max(phi_n - phi_t, 0)) ...
                    * (abs(Delta_t) / D_tMax)^beta_r * Delta_t / abs(Delta_t);

            end
        end

        Tc = [T_t; T_n];
        D = [D_tt, D_tn;
            D_nt, D_nn];

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