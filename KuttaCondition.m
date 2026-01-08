%function to enforce the kutta condition by way of a wake panel which will
%be used to prescribe the doublet strengths at the upper and lower trailing
%edge panels, the result will be influence coefficients mapped to an exact
%solution of the matrix equation

function [L, U_normal] = KuttaCondition(L, S, Xb, Yb, xc, yc, alphaR, c, numPan, n_hat, U_normal)
Sw = 0.000; % Compute wake panel length (currently 3 chord lengths)
Xw1 = Xb(1); % First wake x-coordinate
Xw2 = Xw1 + Sw*cos(alphaR); % Second wake x-coordinate
Yw1 = Yb(1); % First wake y-coordinate
Yw2 = Yw1 + sin(alphaR); % Second wake y-coordinate
xcw = 0.5*(Xw1+Xw2); % Compute wake control point x-coordinate
ycw = 0.5*(Yw1+Yw2); % Compute wake control point y-coordinate

for i = 1:numPan
    dxw = xc(i) - xcw; % Compute x-distance between the ith control point and wake control point
    dyw = yc(i) - ycw; % Compute y-distance between the ith control point and wake control point
    rw2 = dxw^2 + dyw^2; % Compute the squared distance between the ith contorl point and wake control point

    Qiw = (1/(2*pi)) * [-dyw, dxw] / rw2; % Compute wake doublet velocity kernel
    L(i, numPan+1) = Sw * dot(Qiw, n_hat(:,i)); % Compute wake doublet influence coefficient on the ith panel control point and add to doublet matrix
   
end
test = L(:, numPan+1);
 disp(test)
KuttaRow = zeros(1, numPan+1); % Create N+1 row to add to doublet influence coefficient matrix to make system of equations solveable
KuttaRow(1) = 1; % Fixing value of lower TE doublet strength
KuttaRow(numPan) = -1; % Fixing value of upper TE doublet strength
KuttaRow(numPan+1) = -1; % Fixing value of wake strength
L = [L; KuttaRow]
L(:,1) = L(:,1) - L(:,numPan+1); % 
L(:,numPan) = L(:,numPan) + L(:,numPan+1); %
L(:, numPan+1) = []; 
% KuttaRow = zeros(1, numPan+1); % Generate row to append to doublet influence coefficient matrix to
% 
% % Modify Kutta row values to enforce that mu_upper_TE - mu_lower_TE = mu_wake
% KuttaRow(1) = 1; % Lower trailing edge assigned
% KuttaRow(numPan) = -1; % Upper trailing edge assigned
% KuttaRow(numPan+1) = 1; % Wake assigned 
% 
% % Construct modified influence coefficient matrix with Kutta row included
% L(:, numPan+1) = 0; % Modify doublet influence coefficient matrix to include column for wake influence
% D = zeros(numPan, numPan+1); % Initialize updated doublet normal velocity influence coefficient matrix
% D = [L; KuttaRow]; % Append Kutta condition row to influence coefficient matrix
% 
% % Apply Kutta condition constraint to doublet influence coefficient matrix
% D(:,1) = D(:,1) - D(:,numPan+1); 
% D(:,numPan) = D(:,numPan) + D(:,numPan+1);
% 
% % Reduce order of influence coefficient matrix
% D = D(1:numPan, 1:numPan);
% 
% % Create placeholder zero at TE to make system of equations
% % solveable
% U_normal(numPan+1) = 0; % Last value of normal velocity set to make matrix equations solveable




end