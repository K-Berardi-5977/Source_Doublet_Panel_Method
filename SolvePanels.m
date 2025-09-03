function [lambda, Vt, Cp, Nuemann_check] = SolveSourcePanels(I, J, U, beta, numPan, S, rho)

A = zeros(numPan, numPan); %variable to store Integral terms for linear system of equations

for i = 1:numPan
    for j = 1:numPan
        if (i == j)
            A(i,j) = pi; %normal velocity self-influence of the ith panel
        else
            A(i,j) = I(i,j); %assign coefficient value
        end
    end
end

b = zeros(numPan,1); %normal free-stream terms vector (right-hand side of equation)

for n = 1:numPan
    b(n) = -U*2*pi*cos(beta(n)); %compute normal free-stream terms
end

lambda = A\b; %SOLVE LINEAR SYSTEM OF EQUATIONS

lambda_ds = lambda(:).*S(:); %vector contraining the strengths of each panel
Nuemann_check = sum(lambda_ds); %checking that the Neumann boundary condition is satisfied (i.e., no penetration into or out of the surface)

%========== Surface Velocity and Aerodynamic Loads ==========%

Vt = zeros(numPan,1); %Panel tangential velocity vector
Cp = zeros(numPan, 1); %Panel Pressure coefficient vector

for i = 1:numPan
    source_terms = 0;
    for j = 1:numPan
        source_terms = source_terms + (lambda(j)/(2*pi))*(J(i,j)); %contribution of all other jth panels to tangent velocity of ith panel 
    end
    Vt(i) = U*sin(beta(i)) + source_terms; %surface velocity at ith control point
    Cp(i) = 1-(Vt(i)/U)^2; %pressure coefficient evaluated at ith control point
end
end