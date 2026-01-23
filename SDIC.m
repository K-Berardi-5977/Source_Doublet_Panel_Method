%function to compute the influence coefficients of the source and doublet
%panels in both the velocity potential and velocity formulations

function [J, K, L, M, rc2] = SDIC(xc, yc, S, numPan, n_hat, t_hat)

% === Initialize Temporary Variables === %
dxc = zeros(numPan, numPan); % Matrix to store x-distance between control points
dyc = zeros(numPan, numPan); % Matrix to store y-distance between control points
rc2 = zeros(numPan, numPan); % Matrix to store squared distances between control points

for i = 1:numPan
    for j = 1:numPan
        dxc(i,j) = xc(i) - xc(j); % Compute x-distance between jth control point and ith control point
        dyc(i,j) = yc(i) - yc(j); % Compute y-distance between jth control point and ith control point
        rc2(i,j) = dxc(i,j).^2 + dyc(i, j).^2; % Compute squared distance between jth control point and ith control point
    end
end


% === Initialize Influence Coefficient Matrices === %
J = zeros(numPan, numPan); % Doublet tangent influence coefficient
K = zeros(numPan, numPan); % Source tangent velocity influence coefficient
L = zeros(numPan, numPan); % Doublet normal velocity influence coefficient
M = zeros(numPan, numPan); % Source normal velocity influence coefficient

% === Compute Velocity Influence Coefficients ===

for i = 1:numPan
    for j = 1:numPan
        
        if i == j
            J(i,j) = -0.5*S(j) ; % Compute doublet tangential self-influence coefficient
            K(i,j) = 0; % Compute source tangential self-influence coefficient
            M(i,j) = 0.5*S(j); % Compute source normal self-influence coefficient
            L(i,j) = 0; % Compute doublet normal self-influence coefficient
        else
            Vij = (1/(2*pi)) * [dxc(i,j), dyc(i,j)]/(rc2(i,j)); % Compute source velocity kernel
            Qij = (1/(2*pi)) * [dyc(i,j), -dxc(i,j)]/(rc2(i,j)); % Compute doublet velocity kernel
            K(i,j) = S(j) * dot(Vij, t_hat(:,i)); % Compute source tangential influence coefficient of jth panel on the ith panel
            M(i,j) = S(j) * dot(Vij, n_hat(:,i)); % Compute source normal influence coefficient of jth panel on the ith panel

            
            J(i,j) = S(j) * dot(Qij, t_hat(:,i)); % Compute doublet tangent influence coefficient of jth panel on the ith panel
            L(i,j) = S(j) * dot(Qij, n_hat(:,i)); % Compute doublet normal influence coefficient of the ith panel on the jth panel
    end
    end
end

