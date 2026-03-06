% Separate script to iterate over multiple executions of the panel code to
% generate a trend of lift coefficient vs different angles of attack

clc; clear;

alphaD = linspace(0,10, 22) % Generate angles of attack for iteration (in degrees)
U = 1; % Free stream velocity magnitude
c =1;
aoa_data = load('AbbVD_ClDData.dat') % Data from Theory of Wing Sections for NACA 0012 lift coefficient at 5 aoa
Bp = load('foilData.dat'); % Load airfoil grid (boundary) points from data file
cl = zeros(length(alphaD),1); % Initialize lift coefficient vector
cl2 = zeros(length(alphaD),1);

for N = 1:length(alphaD)
result = Dirichilet_ConstantSourceDoublet(Bp, alphaD(N), U, c) % run panel code usin Nth angle of attack
result_archive = Dirichilet_ConstantSourceDoublet_SemiInfWake_Archive(Bp, alphaD(N), U)
cl(N) = result.cl % extract lift coefficients for N angles of attack
cl2(N) = result_archive.cl
end

figure(3); hold on;
plot(alphaD, cl, '-b') % Plot numerical data from code
plot(alphaD, cl2, '*r')
xlabel('X-Coordinate of Airfoil');
ylabel('Coefficient of Lift (C_l)');
plot(aoa_data(:,1), aoa_data(:,2), '^k') % Plot experimental data
title('Coefficient of Lift vs Angle of Attack (\alpha = 5) NACA 0012');
legend('Berardi Discrete Wake', 'Berardi Semi-Infinite Wake','Abbot & von Doenhoff, Re = 6 \times 10^6', Location='northwest');
axis padded
setaxes
hold off