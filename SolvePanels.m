function [sigma, mu, Vt, Cp, CL, Nuemann_check, D, S, D_ds, S_ds, potential] = SolvePanels(xi, yi, Xj, Yj, phi, L, U, beta, N, rho, alpha)
%% ========== READ ME FIRST ==========

% ===== Inputs =====
% xi = x-coordinates of control points
% yi = y-coordinates of control points
% Xj = x-coordinates of panel boundary points
% Yj = y-coordinates of panel boundary points
% S = Panel Lengths
% N = number of panels
% phi = angles between panels and positive x-axis [rad]
% beta = angles between free-stream and outward facing unit normals of panels [rad]
% alpha = angle of attack [rad]
% U = free-stream velocity magnitude
% rho = density of air (1.2 [kg/m^3])

% ===== Outputs =====
% sigma = Nx1 vector of source strengths 
% mu = Nx1 vector of doublet strengths
% Vt = Nx1 vector of panel surface veclocities
% Cp = Nx1 vector of panel pressure coefficients
% CL = scalar total lift coefficient for the airfoil at the given alpha
% Nuemann_check = scalar check of the no-penetration boundary condition
% D = NxN matrix of doublet potential influence coefficients
% S = NxN matrix of source potential influence coefficients
% D_ds = NxN matrix of doublet influence coefficients on  surface velocity
% S_ds = NxN matrix of source influence coefficients for surface velocity

% ====== Code Basic Flow ===== 
% 1) Initialize output variables as needed (decreases processing time)
% 2) Assign source strengths as equal to free-stream normal component (required to solve linear system of equations
% 3) Perform coordinate transform (integral equations for coefficients are much simpler in local panel coordinates)
% 4) Assign source/doublet influenc coefficients based on Katz & Plotkin equations
% 5) Compute Wake Coefficient (temporarily add an N+1 column to the doublet potential influence coefficient matrix)
% 6) Enfource Kutta Condition using
% 7) Solve Linear System of Equations and Check Boundary Conditions
% 8) Compute Surface Velocity and Aerodynamic Loads

%% ========== (1) INITIALIZE VARIABLES ==========
sigma = zeros(N,1); % source strength vector


D = zeros(N, N+1); % doublet potential influence coefficient matrix
S = zeros(N,N); % source potential influence coefficient matrix

D_ds = zeros(N,N); % doublet surface velocity influence coefficient matrix
S_ds = zeros(N,N); % source surface velocity influence coefficient matrix
dl = zeros(N,1)
potential = zeros(N,1); % panel potential vector
% Vt = zeros(N,1)'; % panel surface velocity vector
Cp = zeros(N,1)'; % panel pressure coefficient vector


rotcw = [cos(phi), sin(phi); -sin(phi), cos(phi)]; % clockwise rotation matrix to be dotted with panel global coordinates

%% ========== (2) PRESCRIBE SOURCE STRENGTHS ==========

for i = 1:N
    sigma(i) = U*(cos(alpha)*cos(phi(i)) - sin(alpha)*sin(phi(i))); % source strength at ith panel = free-stream normal velocity at ith panel
end

%% ========== (3, 4, 5) PERFORM COORDINATE TRANSFORM AND COMPUTE INFLUENCE COEFFICIENTS ==========



for i = 1:N
    for j = 1:N
        % ===== Convert to Local Panel Coordinates ===== %
        dx1 = xi(i) - Xj(j); dy1 = yi(i) - Yj(j); % global x and y distance from first panel boundary point to collocation point
        
        
        dx2 = Xj(j+1) - Xj(j); dy2 =  Yj(j+1) - Yj(j); %global x and y components of jth panel
        
        
        X_X1 = dx1*cos(phi(j)) + dy1*sin(phi(j)) %rotate dx1 and dy1 so they are "flat"
        X_X2 = dx2*cos(phi(j)) + dy2*sin(phi(j)) %rotate dx2 and dy2 so they are "flat"

        if i == 1
            dl(j) = X_X2;
        end

        Yo = -dx1*sin(phi(j)) + dy1*cos(phi(j)) -0.01*L(j);
        
        r1 = hypot(X_X1,Yo);
        r2 = hypot((X_X2-X_X1), Yo);

        th1 = atan2(Yo,X_X1);
        th2 = atan2(Yo, (X_X2-X_X1));

        % ===== COMPUTE INFLUENCE COEFFICIENTS =====%
        if i == j 
            D(i,j) = 0.5; %doublet self influence coefficient
            S(i,j) = (1/pi)*(X_X1*log(r1)); %source self influence coefficient

            D_ds(i,j) = 0; % no ith doublet self influence on surface velocity of ith panel
            S_ds(i,j) = (1/(2*pi))*(log(r1/r2)); % influence of ith source on surface velocity of ith panel
        else
            D(i,j) = (-1/(2*pi))*(th2-th1); % doublet influence coefficient of the jth panel doublet on the ith control point
            S(i,j) = (1/(2*pi))*(X_X1*log(r1) - (X_X1-X_X2)*log(r2) + Yo*(th2-th1)); %influence of the jth panel source on the ith control point

            D_ds(i,j) = (-1/(2*pi))*((Yo/(r1^2)) - (Yo/(r2^2)));
            S_ds(i,j) = (1/(2*pi))*(log(r1/r2));
    end
    end
    dxw = xi(i)-Xj(N);
    dyw = yi(i)-Yj(N);
    thw = -atan(dyw/dxw);

    D(i,N+1) = (-1/(2*pi))*thw;
end

%% ========== (6) Enforce Kutta Condition at Trailing Edge ==========

% ==== Add/Substract Wake Doublet Influence Coefficient at Leading/Traling Edges as appliable ===== %
for i = 1:N
    for j = 1:N 
        if j == 1
            D(i,j) = D(i,j) + D(i, N+1);
        elseif j == N
            D(i,j) = D(i,j) - D(i, N+1);
        else
            continue
        end
    end
end

% ===== Remove column N+1 from doublet influence coefficient matrix
D(:,N+1) = [];

%% ========== (7) SOLVE LINEAR SYSTEM OF EQUATIONS AND CHECK BOUNDARY CONDITIONS ==========

% ===== Multiply Source Strenghs by Influence Coefficients and Sum into Right Hand Side ===== %

for i = 1:N 
    for j = 1:N 
        sig_S(i,j) = sigma(j)*S(i,j); % Creates Matrix of source strenghts multiplied by their influence coefficients for each point
    end
end

RHS = sum(sig_S, 2); %right hand Side of Linear System of Equations

% ===== Solve for Doublet Strengths ===== 

mu = (-RHS\D) %solve for doublet strengths

Nuemann_check = sum(sigma(:).*L(:)) %vector contraining the strengths of each panel

%% ========== COMPUTE SURFACE VELOCITY AND AERODYNAMIC LOADS ==========                                                              

% ===== Compute Potential at Each Control Point =====
for i = 1:N
    potential(i) = xi(i)*cos(alpha) + yi(i)*sin(alpha) + dot(mu, D(i,:)) + dot(sigma, S(i,:)); % potential at ith control point
end
% ===== Surface Velocity and Pressure Coefficient =====

for i = 1:N
    s_terms = 0;
    d_terms = 0;
    for j = 1:N
        s_terms = s_terms + sigma(i)*S_ds(i,j);
        d_terms = d_terms + mu(i)*D_ds(i,j);
    end
    Vt(i) = U*(sin(beta(i))) + s_terms + d_terms;
    Cp(i) =  1-(Vt(i)/U)^2;
end

CL = 0

% dphi = diff(mu)
% 
% 
% dphi\dl
end
