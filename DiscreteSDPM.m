 %========= NUMERICAL SOURCE PANEL METHOD CODE =========%
clc; 
clear;

%========== Define Knowns ==========%
U = 1; %free stream velocity
alphaD =5; %angle of attack in degrees
alphaR = alphaD*(pi/180); %angle of attack in radians
c = 1; %chord length of foil;
rho = 1.2; %density of air (kg/cubic meter)

%========== Define Airfoil Profile ==========%
def_foil = 'Use .dat File'; %variable to control how profile is generated

[XB, YB, XC, YC, phiR, betaR, S, numPan] = LoadPanels(def_foil, c, alphaD);


%======= Determine Normal & Tangentential Infuence Coefficients =======%
% [I, J] = SDPM_InfluenceCoeff(XC, YC, XB, YB, phiR, S, numPan); %function solving for influence coefficients

%========== Solve Linear System of Equations ==========%
[lambda, mu, Vt, Cp, Neumann_check, A, B, potential, RHS, DL] = SolvePanels(XC, YC, XB, YB, S, phiR, U, betaR, numPan, rho, alphaR);

%========== Plot Streamlines ==========%
% [Nxx , Nyy, Vxy, theta_pj, psi, THETA, Cpxy_mask] = PM_streamlines(XC, YC, XB, YB, phiR, S, lambda, U, alphaD, Cp, numPan);

figure; hold on; axis equal;
plot(XB,YB, 'b.', MarkerSize=7);
plot(XC, YC, 'r*');
plot(XB,YB,'k');
% plot(x_c(indices), y_c(indices), 'bo', MarkerSize=7, MarkerFaceColor='c')
title('Discretized Body Panels')
xlabel('X')
ylabel('Y')
legend('Panel Bounds', 'Control Points')



XB(end) = [];
half_x = floor(numPan/2);
figure; hold on;
set(gca, 'YDir','reverse')
plot(XB(1:half_x), Cp(1:half_x), 'b');
plot(XB(half_x+1:end), Cp(half_x+1:end)); 
% plot(x_c, V_s, 'r--');
title(['Pressure Distribution on Airfoil Surface ($\alpha = ', num2str(alphaD), ')$'], 'Interpreter','latex');
xlabel('X-Coordinate of Airfoil');
ylabel('Coefficient of Pressure (Cp)');
legend('Bottom Cp', 'Top Cp');




