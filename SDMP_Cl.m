% Separate script to iterate over multiple executions of the panel code to
% generate a trend of lift coefficient vs different angles of attack

clc; clear;

aoa_data = load('AbbVD_ClDData.dat'); % Data from Theory of Wing Sections for NACA 0012 lift coefficient at 5 aoa
% aoa_data = load('xfoil_alpha_cl.dat');
alphaD = aoa_data(:,1); % Generate angles of attack for iteration (in degrees)
U = 1; % Free stream velocity magnitude
c =1;
rho = 1.225;

cm_data = load('xfoil_alpha_cm_4412.dat')
[Bp, Nw] = call_naca4; % Load airfoil grid (boundary) points from data file
cl = zeros(length(alphaD),1); % Initialize lift coefficient vector
cm = zeros(length(alphaD),1);
cl2 = zeros(length(alphaD),1);

for N = 1:length(alphaD)
result = Dirichilet_ConstantSourceDoublet(Bp, Nw, alphaD(N), U, c, rho) % run panel code usin Nth angle of attack
result_archive = Dirichilet_ConstantSourceDoublet_SemiInfWake_Archive(Bp, alphaD(N), U)
cl(N) = result.cl % extract lift coefficients for N angles of attack
cm(N) = result.aero.cm;
cl2(N) = result_archive.cl
end

figure(3); hold on;
plot(alphaD, cl, 'ob') % Plot numerical data from code
% plot(alphaD, cl2, '*r') % Plot numerical data from code
% xlabel('\alpha');
% ylabel('Coefficient of Lift (c_l)');
plot(aoa_data(:,1), aoa_data(:,2), '^k') % Plot experimental data

% title('Coefficient of Lift vs Angle of Attack for NACA 0012 Airfoil');
legend('Berardi SDPM Code', 'Abbott & von Doenhoff, Re = 6 \times 10^6', Location='northwest');
axis padded
plot(alphaD, cl, '-b', LineWidth=0.75, HandleVisibility='off') % Plot numerical data from code
setaxes
hold off

% figure(4); hold on;
% plot(alphaD, cl, 'ob') % Plot numerical data from code
% xlabel('\alpha');
% ylabel('Coefficient of Lift (c_l)');
% plot(aoa_data(:,1), aoa_data(:,2), '*r') % Plot experimental data
% title('Coefficient of Lift vs Angle of Attack for NACA 4412 Airfoil');
% legend('Berardi Discrete Wake','Xfoil Results, Re = 5 \times 10^4', Location='northwest');
% axis padded
% setaxes
% hold off
% figure(4); hold on;
% plot(alphaD, cm, '-b') % Plot numerical data from code
% plot(cm_data(:,1),cm_data(:,2), '^k');
% ylabel('Coefficient of Moment (c_m)');
% title('Moment vs Angle of Attack ');
% legend('Berardi Discrete Wake', 'Xfoil Results, Re = 5 \times 10^4', Location='northwest');
% axis padded
% setaxes
% hold off