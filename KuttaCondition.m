%function to enforce the kutta condition by way of a wake panel which will
%be used to prescribe the doublet strengths at the upper and lower trailing
%edge panels, the result will be influence coefficients mapped to an exact
%solution of the matrix equation

function [D, U_normal] = KuttaCondition(L, numPan, U_normal)

KuttaRow = zeros(1, numPan+1); % Generate row to append to doublet influence coefficient matrix to

% Modify Kutta row values to enforce that mu_upper_TE - mu_lower_TE = mu_wake
KuttaRow(1) = 1; % Lower trailing edge assigned
KuttaRow(numPan) = -1; % Upper trailing edge assigned
KuttaRow(numPan+1) = 1; % Wake assigned 

% Construct modified influence coefficient matrix with Kutta row included
L(:, numPan+1) = 0; % Modify doublet influence coefficient matrix to include column for wake influence
D = zeros(numPan, numPan+1); % Initialize updated doublet normal velocity influence coefficient matrix
D = [L; KuttaRow]; % Append Kutta condition row to influence coefficient matrix

% Apply Kutta condition constraint to doublet influence coefficient matrix
D(:,1) = D(:,1) - D(:,numPan+1); 
D(:,numPan) = D(:,numPan) + D(:,numPan+1);

% Reduce order of influence coefficient matrix
D = D(1:numPan, 1:numPan);

% Create placeholder zero at TE to make system of equations
% solveable
U_normal(numPan+1) = 0; % Last value of normal velocity set to make matrix equations solveable




end