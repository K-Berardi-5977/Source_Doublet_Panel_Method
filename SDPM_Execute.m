% ===== SOURCE-DOUBLET PANEL CODE TEST SCRIPT ===== %
% this code is for running the function and is primarily intended to create
% Cp comparisons -- Cl comparisons are generated with SDMP_Cl.m 

% In its current state many of the data comparisons are hard coded in and
% are only valid when the correct angle of attack is selected based on the
% comment label 

clc
clear
%% ========== INPUT PARAMETERS ========= %
alphaD = 5; % Angle of attack (degrees)
U = 1; % Free stream velocity magnitude
c = 1; % Chord Length (m)
m = 1.48E-5; % kinematic viscosity of air at 15 degrees C
rho = 1.2; % Air density (kg/ cubic meter)

% Cp vs x/c for 10 degree angle of attack
ExpData1 = load('Cp_Gregory_Oreilly.dat');
ExpData2 = load('Cp_Ladson.dat')

% Cp vs x/c for 5 degree angle of attack
ExpData3 = load('KPData.dat') % Katz and Plotkin Code -- Cp vs x/c @ aoa = 5

%% ===== Initialize GEOMETRY ===== %
Bp = load('foilData.dat'); % Load airfoil grid (boundary) points from data file

% ===== Perform Computations ===== %
result = Dirichilet_ConstantSourceDoublet(Bp, alphaD, U)

% ===== Plot Results =====
% READ LABELS ABOVE FIGURE PLOTS

% For alphaD = 10 use this portion and comment out ExpData3 and figure(2)
% figure(1); hold on;
% set(gca, 'YDir','reverse')
% plot(result.X_Cp, result.Cp, '-b')
% plot(ExpData1(:,1), ExpData1(:,2), 'r*')
% plot(ExpData2(:,1), ExpData2(:,2), 'go')
% title(['Pressure Distribution on Airfoil Surface ($\alpha = ', num2str(alphaD), ')$'], 'Interpreter','latex');
% xlabel('X-Coordinate of Airfoil');
% ylabel('Coefficient of Pressure (Cp)');
% legend('Berardi SDPM Code', 'Gregory, Re = 3 \times 10^7', 'Ladson, Re = 6 /times 10^7', 'Katz & Plotkin SDPM Code', Location='northeast')
% axis padded
% hold off

% For alphaD = 5 use this portion and comment out ExpData1, ExpData2, and figure(1)
figure(2); hold on;
set(gca, 'YDir','reverse')
plot(result.X_Cp, result.Cp, '-b')
plot(ExpData3(:,1), ExpData3(:,2), '^k')
title(['Pressure Distribution on Airfoil Surface ($\alpha = ', num2str(alphaD), ')$'], 'Interpreter','latex');
xlabel('X-Coordinate of Airfoil');
ylabel('Coefficient of Pressure (Cp)');
legend('Berardi SDPM Code', 'Katz and Plotkin SDPM Code', Location='northeast')
axis padded
hold off