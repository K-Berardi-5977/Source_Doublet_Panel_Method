% Function to prescribe source strength per the no internal perturbation Dirichlet Boundary Condition

% === Inputs ===
% Bp - panel boundary points as (x,y) pairs ordered clockwise
% alphaR - pitch angle (angle of attack) [radians]
% cth = mesh.cos_th - cos(th) for each panel (useful for panel frame rotations) (numPan x 1)
% sth = mesh.sin_th - sin(th) for each panel (useful for panel frame rotations) (numPan x 1)

% === Outputs ===
% sigma - body panel source strengths (numPan x 1)

function sigma = PreSource(U, alphaR, mesh)

    sigma = U*(cos(alphaR).*mesh.sin_th - sin(alphaR).*mesh.cos_th); % Constant-strength sources equal to freestream normal component on each panel.

end