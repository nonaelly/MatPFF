% PPR cohesive model
Delta_t=1
Delta_n=1
delta_n=1
delta_t=1
psi = min(phi_n, phi_t) +

%% Recognition process
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


    if abs(Delta_t) >= d_tMax
        % soften
        D_nn = Gamma_n / delta_n^2 * ((m^2 - m) * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^(m - 2) ...
            + (alpha_k^2 - alpha_k) * (1 - Delta_n / delta_n)^(alpha_k - 2) * (m / alpha_k + Delta_n / delta_n)^m ...
            - 2 * alpha_k * m * (1 - Delta_n / delta_n)^(alpha_k - 1) * (m / alpha_k + Delta_n / delta_n)^(m - 1))...
            * (Gamma_t * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^n + max(phi_t - phi_n, 0));

        D_nt = Gamma_n * Gamma_t / (delta_n * delta_t) * (m * (1 - Delta_n / delta_n)^alpha_k * (m / alpha_k + Delta_n / delta_n)^(m - 1) ...
            - alpha_k * (1 - Delta_n / delta_n)^(alpha_k - 1) * (m / alpha_k + Delta_n / delta_n)^m)...
            * (n * (1 - abs(Delta_t) / delta_t)^beta_k * (n / beta_k + abs(Delta_t) / delta_t)^(n - 1) ...
            - beta_k * (1 - abs(Delta_t) / delta_t)^(beta_k - 1) * (n / beta_k + abs(Delta_t) / delta_t)^n) * Delta_t / abs(Delta_t);

        T_t = Gamma_t / delta_t * (n * (1 - Delta_t / delta_t)^beta_k * (n / beta_k + Delta_t / delta_t)^(n - 1) ...
            - beta_k * (1 - Delta_t / delta_t)^(beta_k - 1) * (n / beta_k + Delta_t / delta_t)^n)...
            * (Gamma_n * (1 - abs(Delta_n) / delta_n)^alpha_k * (m / alpha_k + abs(Delta_n) / delta_n)^m + max(phi_n - phi_t, 0)) * Delta_t / abs(Delta_t);

    elseif abs(Delta_t) < d_tMax
        % unload / reload

    end
end