%function to enforce the kutta condition by way of a wake panel which will
%be used to prescribe the doublet strengths at the upper and lower trailing
%edge panels, the result will be influence coefficients mapped to an exact
%solution of the matrix equation

function [D] = KuttaCondition(L, numPan)

KuttaRow = zeros(1, numPan+1); % Generate row to append to doublet influence coefficient matrix to

% Modify Kutta row values to enforce that mu_upper_TE - mu_lower_TE = mu_wake
KuttaRow(1) = 1; % Lower trailing edge assigned
KuttaRow(numPan) = -1; % Upper trailing edge assigned
KuttaRow(numPan+1) = 1; % Wake assigned 

% Construct modified influence coefficient matrix with Kutta row included
D = zeros(numPan, numPan+1); % Initialize updated doublet normal velocity influence coefficient matrix
D = [D; KuttaRow]; % Append Kutta condition row to influence coefficient matrix



end