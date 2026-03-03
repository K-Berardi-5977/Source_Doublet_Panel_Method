% Function to preprocess inputs:

% === INPUTS ===
% Bp - panel boundary points as (x,y) pairs ordered counter clockwise
% alphdaD - pitch angle (angle of attack) [degrees]

% === Outputs ===
% Bp - panel boundary points as (x,y) pairs ordered clockwise
% NN - number of boundary points
% numPan - number of body panels
% num_d - number of doublets (body + wake)
% alphaR - pitch angle (angle of attack) [radians]
% diag_idx - diagonal index variable for any numPan x numPan matrix

function [Bp, alphaR, NN, numPan, num_d, diag_idx] = preprocess(Bp, alphaD)

Bp = flipud(Bp);                 % Force clockwise paneling
alphaR = (pi/180)*alphaD;        % Convert pitch angle to radians

NN     = size(Bp, 1);            % Number of panel boundary points
numPan = NN - 1;                 % Nnumber of body panels
num_d  = NN;                     % Number of wake panels

% Diagonal indexing for any (numPan x numPan) matrix (primarily for
% self-influence coefficients)
diag_idx = 1:(numPan+1):(numPan^2);

end