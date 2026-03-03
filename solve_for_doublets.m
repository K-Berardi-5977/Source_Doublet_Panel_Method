% === INPUTS ===
% sys - struct containing matrix that entirely defines the zero internal perturbation potential Dirichlet condition matrix equations 
% for potential, and contains wake panel coordinates 


% === OUTPUTS ===
% mu - vector of doublet strengths (including wake)

function mu = solve_for_doublets(sys)
    b  = sys.A_full(:, end); % RHS column
    A  = sys.A_full(:, 1:end-1); % doublet coefficients (body + wake)
    mu = A \ b; % Solve for doublet strengths (validation: mu(end-1) - mu(1) = mu(end) if dont correctly)
end
