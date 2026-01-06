% ===== SOURCE-DOUBLET PANEL CODE REV 2 ===== %
clc
clear
%% ========== INPUT PARAMETERS ========= %
alphaD = 10; % Angle of attack (degrees)
alphaR = alphaD*(pi/180); % Angle of attack (radians)
U = 1; % Free stream velocity magnitude (m/s)
U_inf = U*[cos(alphaR); sin(alphaR)] % Free stream velocity vector (m/s)
c = 1; % Chord Length (m) 
rho = 1.2; % Air density (kg/ cubic meter)

%% ========== PRE-PROCESSING ========== %

% ===== Load Airfoil Geometry (Panel Mesh) ===== %
[XB, YB, XC, YC, phiR, betaR, S, numPan, n, t] = LoadPanels2(alphaD);