function [Bp, Nw] = call_naca4()
%PROMPT_NACA4_AIRFOIL Prompt user for airfoil parameters and generate Bp
%
% This wrapper asks the user for NACA airfoil parameters and calls
% naca4_airfoil() to generate the boundary points.

fprintf('\n=== NACA 4-Digit Airfoil Generator ===\n');

% ---- Airfoil designation ----
digits = input('Enter NACA 4-digit code (e.g. 0012, 2412): ','s');

% ---- Chord length ----
chord = input('Enter chord length: ');

% ---- Panel count ----
numPan = input('Enter total number of panels (must be even): ');

if mod(numPan,2) ~= 0
    error('Number of panels must be even.');
end

% ---- Spacing type ----
spacingType = input('Spacing type ("cosine" or "uniform"): ','s');

% ---- Trailing edge option ----
teOption = input('Closed trailing edge? (1 = yes, 0 = no): ');

closedTE = logical(teOption);

% ---- Wake panels ----
Nw = round(input('Enter starting number of wake panels: '));

if Nw < 1
    error('Number of wake panels must be >= 1.');
end

% ---- Generate boundary points ----
Bp = naca4_airfoil(digits, chord, numPan, spacingType, closedTE);

fprintf('\nAirfoil successfully generated.\n');
fprintf('Boundary points: %d\n', size(Bp,1));
fprintf('Panels: %d\n\n', numPan);

end