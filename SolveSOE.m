% Function to solve the system of matrix equations for the source doublet
% panel method and return a vector of doublet strengths and panel surface
% velocities

function [mu] = SolveSOE(L, M, U_normal, sigma, numPan)

S = M * sigma'; % Compute source terms by multiplying source strengths by their influence coefficients 

for i = 1:numPan
    % Create right hand side of matrix equation by subtracting the free-stream and source terms
            % if i == numPan
            RHS(i,1) = - S(i);
           
end
% === Solve System of Equations for Doublet Strength === %
mu = L\RHS;
 
end