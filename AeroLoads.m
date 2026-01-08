function [Cp, VT, Vt_s, Vt_d] = AeroLoads(U, U_tangent, S, sigma, mu, J, K, numPan)

% === Compute Singularity Tangent Velocity Terms and Pressure Coefficient === %
Vt_s = sigma*K; % Compute source tangent velocity contribution at each panel
Vt_d = mu*J; % Compute doublet tangent velocity contribution at each panel
VT = zeros(numPan, 1); % Create vector to store panel surface velocities
Cp = zeros(numPan, 1); % Create vector to store pressure coefficient vallues

for i = 1:numPan
    VT(i) = U_tangent(i) + Vt_s(i) + Vt_d(i); % Compute total tangent velocity at ith panel
    Cp(i) = 1-(VT(i)/U)^2; % Compute pressure coefficient of ith panel
end
Gamma = 0;
end