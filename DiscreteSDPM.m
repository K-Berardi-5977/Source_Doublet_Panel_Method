 %========= NUMERICAL SOURCE PANEL METHOD CODE =========%
clc; 
clear;

%========== Define Knowns ==========%
U = 1; %free stream velocity
alphaD =10; %angle of attack in degrees
alphaR = alphaD*(pi/180); %angle of attack in radians
c = 1; %chord length of foil;
rho = 1.2; %density of air (kg/cubic meter)

%========== Define Airfoil Profile ==========%
def_foil = 'Use .dat File'; %variable to control how profile is generated

[XB, YB, XC, YC, phiR, betaR, L, numPan] = LoadPanels(def_foil, c, alphaD);


%======= Determine Normal & Tangentential Infuence Coefficients =======%
% [I, J] = SDPM_InfluenceCoeff(XC, YC, XB, YB, phiR, S, numPan); %function solving for influence coefficients

%========== Solve Linear System of Equations ==========%
[sigma, mu, Vt, Cp, CL, Nuemann_check, D, S, D_ds, S_ds, potential] = SolvePanels(XC, YC, XB, YB, phiR, L, U, betaR, numPan, rho, alphaR);

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
plot(XC(1:half_x), Cp(1:half_x), 'b');
plot(XC(1:half_x), Cp(1:half_x), 'bo')
plot(XC(half_x:end), Cp(half_x:end), 'r');
plot(XC(half_x:end), Cp(half_x:end), 'ro');
% plot(x_c, V_s, 'r--');
title(['Pressure Distribution on Airfoil Surface ($\alpha = ', num2str(alphaD), ')$'], 'Interpreter','latex');
xlabel('X-Coordinate of Airfoil');
ylabel('Coefficient of Pressure (Cp)');
legend('Bottom Cp', 'Top Cp');




