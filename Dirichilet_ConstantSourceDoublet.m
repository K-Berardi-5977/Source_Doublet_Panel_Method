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

function result = Dirichilet_ConstantSourceDoublet(Bp, alphaD, U, c)

 %% ===== 1) Preprocess Data =====
    [Bp, alphaR, NN, numPan, num_d, diag_idx, Nw, Lw] = preprocess(Bp, alphaD, c);

    %% ===== 2) Construct Body Mesh =====
    mesh = build_panels(Bp, numPan);
    
    %% ===== 3) Construct Wake Mesh =====
    wake = build_wake(mesh, U, alphaR, c, Nw, Lw);

    %% ===== 4) Fix the source strengths on the surface =====
    sigma = PreSource(U, alphaR, mesh);

    %% ===== 5) Compute Influence Coefficient Matrices =====
    inflBB = SDPM_Influence(mesh, wake, numPan, diag_idx); % Influence coefficients for body panels
    inflBW = SDPM_Influence_BodyWake(mesh, wake); % Influence coefficients of wake panels

    %% ===== 6) Assemble System of Equations and Solve for Doublet Panel Strengths =====
    sys  = Internal_Dirichlet_Formulation(inflBB, inflBW, sigma, mesh, numPan, num_d);
    mu   = solve_for_doublets(sys);

    %% ===== 7) Postprocess Steps 2-5 to Compute Surface Velocity and Aerodynamic Load Coefficients =====
    aero = postprocess_aero(mu, sigma, inflBB, inflBW, mesh, U, alphaR, numPan, num_d, sys);

    %% ===== 9) Wrap results in Struct =====
    result = ExportResults(mu, aero, mesh, numPan);
  
    
   



end

