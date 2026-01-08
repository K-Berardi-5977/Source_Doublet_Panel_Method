% ===== SOURCE-DOUBLET PANEL CODE REV 2 ===== %
clc
clear
%% ========== INPUT PARAMETERS ========= %
alphaD = 9; % Angle of attack (degrees)
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

[J, K, L, M, rc2] = SDIC(xc, yc, S, numPan, n_hat, t_hat);

% ===== Enforce Kutta Condition ===== %

[L, U_normal] = KuttaCondition(L, S, Xb, Yb, xc, yc, alphaR, c, numPan, n_hat, U_normal);

% ===== Solve System of Equations ===== %

[mu] = SolveSOE(L, M, U_normal, sigma, numPan)

[Cp, VT, Vt_s, Vt_d] = AeroLoads(U, U_tangent, S, sigma, mu, J, K, numPan);

Xb(end) = [];
half_x = floor(numPan/2);
figure; hold on;
set(gca, 'YDir','reverse')
plot(xc(1:half_x), Cp(1:half_x), 'b');
plot(xc(half_x:end), Cp(half_x:end), 'r');
plot(xc(1:half_x), Cp(1:half_x), 'bo')
plot(xc(half_x:end), Cp(half_x:end), 'ro');
% plot(xc, VT)
% plot(x_c, V_s, 'r--');
title(['Pressure Distribution on Airfoil Surface ($\alpha = ', num2str(alphaD), ')$'], 'Interpreter','latex');
xlabel('X-Coordinate of Airfoil');
ylabel('Coefficient of Pressure (Cp)');
legend('Bottom Cp', 'Top Cp');