% ===== SOURCE-DOUBLET PANEL CODE TEST SCRIPT ===== %
% this code is for running the function and is primarily intended to create
% Cp comparisons -- Cl comparisons are generated with SDMP_Cl.m 

% In its current state many of the data comparisons are hard coded in and
% are only valid when the correct angle of attack is selected based on the
% comment label 

clc
clear
%% ========== INPUT PARAMETERS ========= %
alphaD = 10; % Angle of attack (degrees)
U = 1; % Free stream velocity magnitude
c = 1; % Chord Length (m)
rho = 1.225; % Air density (kg/ cubic meter)

% Cp vs x/c for 10 degree angle of attack
ExpData1 = load('Cp_Gregory_Oreilly.dat');
ExpData2 = load('Cp_Ladson.dat');

% Cp vs x/c for 5 degree angle of attack
ExpData3 = load('KPData.dat'); % Katz and Plotkin Code -- Cp vs x/c @ aoa = 5

%% ===== Initialize GEOMETRY ===== %
Bp = load('foilData.dat'); % Load airfoil grid (boundary) points from data file
[Bp, Nw] = call_naca4;

% ===== Perform Computations ===== %
result = Dirichilet_ConstantSourceDoublet(Bp, Nw, alphaD, U, c, rho);
result_archive = Dirichilet_ConstantSourceDoublet_SemiInfWake_Archive(Bp, alphaD, U)

% ===== Plot Results =====
% READ LABELS ABOVE FIGURE PLOTS

% For alphaD = 10 use this portion and comment out ExpData3 and figure(2)
figure(1); hold on;
set(gca, 'YDir','reverse')
plot(result.X_Cp, result.Cp, 'ob')
plot(result.X_Cp, result.Cp, '-b', LineWidth=0.85, Handlevisibility = 'off')
% plot(result_archive.X_Cp, result_archive.Cp,'*r')
plot(ExpData1(:,1), ExpData1(:,2), '^k', LineWidth=1.2)
plot(ExpData2(:,1), ExpData2(:,2), 'r*')
% title(['Pressure Distribution on Airfoil Surface (\alpha = ', num2str(alphaD), ')']);
% xlabel('x/c');
% ylabel('Coefficient of Pressure (Cp)');
legend('Berardi SDPM Code','Gregory, Re = 3 \times 10^7', 'Ladson, Re = 6 \times 10^7', Location='northeast')
setaxes
axis padded
hold off

% For alphaD = 5 use this portion and comment out ExpData1, ExpData2, and figure(1)
% figure(2); hold on;
% set(gca, 'YDir','reverse')
% plot(result.X_Cp, result.Cp, '-b')
% % plot(result_archive.X_Cp, result_archive.Cp, '*r')
% plot(ExpData3(:,1), ExpData3(:,2), '^k')
% title(['Pressure Distribution on Airfoil Surface (\alpha = ', num2str(alphaD), ')']);
% xlabel('X-Coordinate of Airfoil');
% ylabel('Coefficient of Pressure (Cp)');
% legend('Beardi Discrete Wake', 'Katz and Plotkin Numerical Diff', Location='northeast')
% axis padded
% setaxes
% hold off

% figure(3); hold on;
% set(gca, 'YDir','reverse')
% plot(Bp, result.Cp, '-b')
% % plot(result_archive.X_Cp, result_archive.Cp, '*r')
% plot(ExpData3(:,1), ExpData3(:,2), '^k')
% title(['Pressure Distribution on Airfoil Surface (\alpha = ', num2str(alphaD), ')']);
% xlabel('X-Coordinate of Airfoil');
% ylabel('Coefficient of Pressure (Cp)');
% legend('Beardi Discrete Wake', 'Katz and Plotkin Numerical Diff', Location='northeast')
% axis padded
% setaxes
% hold off

