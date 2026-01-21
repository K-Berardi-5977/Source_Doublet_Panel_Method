function [Cp, VT, Vt_s, Vt_d, Vt_w] = AeroLoads(U, U_tangent, S, sigma, mu, J, Jwake, K, numPan)

% === Compute Singularity Tangent Velocity Terms and Pressure Coefficient === %
Vt_s = K * sigma; % Compute source tangent velocity contribution at each panel
Vt_d = J * mu; % Compute doublet tangent velocity contribution at each panel
Vt_w = Jwake * (mu(1)-mu(numPan));
VT = zeros(numPan, 1); % Create vector to store panel surface velocities
Cp = zeros(numPan, 1); % Create vector to store pressure coefficient vallues

VT = U_tangent + Vt_s + Vt_d + Vt_w; % Compute total tangent velocity at ith panel
Cp = 1-((VT)/U).^2;
end