% This function provides the lift and pressure distribution characteristics
% for a 2-D airfoil in an inviscid, incompressible, irrotational flow

% 2-D Dirichlet panel method: constant-strength sources + constant-strength doublets
% with a single wake doublet (Kutta enforced at trailing edge (TE))

%{INPUTS}%
% Bp = (X,Y) Cartesion boundary points of body panels
% alphaD = angle of attack in degrees
% U = free-stream velocity magnitude

%{OUTPUTS}%
% result.mu = vector containin the doublet strength of each panel
% result.Cp = pressure coefficient at the control point of each panel
% result.cl = lift coefficient based on circulation and Kelvin's theorem of lift
% result.X_Cp = control point x-coordinates for pressure coefficient plotting

function result = Dirichilet_ConstantSourceDoublet(Bp, alphaD, U)

 %% ===== 1) Preprocess Data =====
    [Bp, alphaR, NN, numPan, num_d, diag_idx] = preprocess(Bp, alphaD);

    %% ===== 2) Construct Mesh =====
    mesh = build_panels(Bp, numPan);

    %% ===== 3) Fix the source strengths on the surface =====
    sigma = PreSource(U, alphaR, mesh);

    %% ===== 4) Compute Influence Coefficient Matrices =====
    infl = SDPM_Influence(mesh, numPan, diag_idx);

    %% ===== 5) Assemble System of Equations and Solve for Doublet Panel Strengths =====
    sys  = Internal_Dirichlet_Formulation(infl, sigma, mesh, numPan, num_d);
    mu   = solve_for_doublets(sys);

    %% ===== 6) Postprocess Steps 2-5 to Compute Surface Velocity and Aerodynamic Load Coefficients =====
    aero = postprocess_aero(mu, sigma, infl, mesh, U, alphaR, numPan, num_d, sys);

    %% ===== 7) Wrap results in Struct =====
    result = ExportResults(mu, aero, mesh, numPan);


% % Initialize variables
% phi = cx*cos(alphaR) + cz*sin(alphaR) + mu(1:numPan); % Compute panel potentials
% Cp2 = zeros(numPan-1, 1); % Initialize Pressure Coefficient Matrix
% Vt = zeros(numPan-1, 1); % Initialize panel surface velocity vector
% 
% % Compute surface velocity using forward difference method
% dS = (Lp(2:numPan) + Lp(1:numPan-1))/2; % Surface increment
% dphi = phi(1:numPan-1) - phi(2:numPan); % Potential increment
% Vt2 = dphi ./ dS; % Surface velocity (surface derivative of potential)
% Cp2 = 1 - (Vt2./U).^2; % Pressure coefficient (steady, inviscid, incompressible)
% 


end

