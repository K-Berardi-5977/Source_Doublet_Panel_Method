% ===== SOURCE-DOUBLET PANEL CODE REV 2 ===== %
clc
clear
%% ========== INPUT PARAMETERS ========= %
alphaD = 10; % Angle of attack (degrees)
alphaR = alphaD*(pi/180); % Angle of attack (radians)
U = 1; % Free stream velocity magnitude (m/s)
U_inf = U*[cos(alphaR); sin(alphaR)]; % Free stream velocity vector (m/s)
c = 1; % Chord Length (m) 
rho = 1.2; % Air density (kg/ cubic meter)

%% ========== PRE-PROCESSING ========== %

% ===== Load Airfoil Geometry (Panel Mesh) ===== %
[Xb, Yb, xc, yc, betaR, S, numPan, n_hat, t_hat] = LoadPanels2(alphaD);

% ===== Rotate Free-Stream into Local Panel Coordinate Frame ===== %
[U_tangent, U_normal, U_local] = LocalFreeStream(U_inf, numPan, t_hat, n_hat);




%% ========== SIMULATION EXECUTION ========== %

% ===== Initialize Variables and Enforce Fixed Parameters ===== %

[sigma, Neumann_Validation] = FixSourceStrength(U_normal, S);

[J, K, L, M, rc2] = SDIC(xc, yc, S, numPan, n_hat);

% ===== Enforce Kutta Condition ===== %

[D, U_normal] = KuttaCondition(L, numPan, U_normal)

% ===== Solve System of Equations ===== %

[mu] = SolveSOE(D, M, U_normal, sigma, numPan)