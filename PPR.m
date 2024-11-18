% PPR cohesive model
sigma_max = 0.5; % MPa
m = 0.1672;
alpha_k = 3.10;
delta_n = 2.5; % mm
tau_max = 0.5; % MPa
n = 0.1672;
beta_k = 3.10;
delta_t = 2.5; % mm

%% Recognition process
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
    * (Gamma_n * (n / beta_k) ^ n + max(phi_t - phi_n, 0));
% unload / reload parameter
alpha_r = 1;
beta_r = 1;


if Delta_n > delta_n || (Delta_n <= delta_n && abs(Delta_t) > delta_t)
    % interface failure
    D_nn = 1e-15;

    D_nt = 1e-15;

    T_n = 1e-15;
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
        % unload / reload
        D_tt = Gamma_t / delta_t * (n * (1 - D_tMax / delta_t)^beta_k * (n / beta_k + D_tMax / delta_t)^(n - 1) ...
            - beta_k * (1 - D_tMax / delta_t)^(beta_k - 1) * (n / beta_k + D_tMax / delta_t)^n)...
            * (Gamma_n * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^m + max(phi_n - phi_t, 0))...
            * beta_r / D_tMax * (Delta_t / D_tMax)^(beta_r - 1);

        D_tn = Gamma_n * Gamma_t / (delta_n * delta_t) * (m * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^(m - 1) ...
            - alpha_k * (1 - Delta_n / delta_n)^(alpha_k - 1) * (m / alpha_k + Delta_n / delta_n)^m)...
            * (n * (1 - D_tMax / delta_t)^beta_k * (n / beta_k + D_tMax / delta_t)^(n - 1) ...
            - beta_k * (1 - D_tMax / delta_t)^(beta_k - 1) * (n / beta_k + D_tMax / delta_t)^n) * Delta_t / abs(Delta_t) * (Delta_t / D_tMax)^beta_r;

        T_t = Gamma_t / delta_t * (n * (1 - D_tMax / delta_t)^beta_k * (n / beta_k + D_tMax / delta_t)^(n - 1) ...
            - beta_k * (1 - D_tMax / delta_t)^(beta_k - 1) * (n / beta_k + D_tMax / delta_t)^n)...
            * (Gamma_n * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^m + max(phi_n - phi_t, 0)) ...
            * (abs(Delta_t) / D_tMax)^alpha_r;

    end
end