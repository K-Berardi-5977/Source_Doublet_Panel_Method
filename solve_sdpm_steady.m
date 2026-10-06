function result = solve_sdpm_steady(cfg)
%SOLVE_SDPM_STEADY
% Top-level steady source-doublet panel-method solver.
%
% INPUT
%   cfg : configuration structure containing:
%
%       cfg.airfoil.code
%       cfg.airfoil.chord
%       cfg.airfoil.nPanels
%       cfg.airfoil.spacing
%       cfg.airfoil.closedTE
%
%       cfg.flow.U
%       cfg.flow.alphaDeg
%       cfg.flow.rho
%
%       cfg.wake.nPanels
%       cfg.wake.lengthChord
%
%       cfg.reference.xMomentChord
%       cfg.reference.zMomentChord
%
% OUTPUT
%   result : structure containing flow, geometry, singularity strengths,
%            aerodynamic outputs, and linear-system information.
%
% SOLVER SEQUENCE
%
%   1) Define derived flow quantities
%   2) Generate NACA geometry
%   3) Build body panel geometry
%   4) Build prescribed steady wake
%   5) Prescribe source strengths
%   6) Compute body-body influence coefficients
%   7) Compute wake-body influence coefficients
%   8) Assemble Dirichlet system
%   9) Solve for doublet strengths
%  10) Postprocess aerodynamic quantities
%  11) Package results


%% ========================================================================
%  1) DEFINE FLOW
% =========================================================================

flow = cfg.flow;

% Convert user input from degrees to radians
flow.alphaRad = deg2rad(flow.alphaDeg);

% Global freestream velocity vector
flow.Uvec = flow.U .* ...
    [cos(flow.alphaRad), sin(flow.alphaRad)];


%% ========================================================================
%  2) GENERATE AIRFOIL GEOMETRY
% =========================================================================

Bp = naca4_airfoil( ...
    cfg.airfoil.code, ...
    cfg.airfoil.chord, ...
    cfg.airfoil.nPanels, ...
    cfg.airfoil.spacing, ...
    cfg.airfoil.closedTE);


%% ========================================================================
%  3) BUILD BODY PANEL GEOMETRY
% =========================================================================

% Generic geometry routine used for any set of connected panels.
body = panel_geometry(Bp);


%% ========================================================================
%  4) BUILD PRESCRIBED STEADY WAKE
% =========================================================================

% Convert nondimensional user wake length to dimensional length.
wakeCfg = cfg.wake;

wakeCfg.length = ...
    wakeCfg.lengthChord .* cfg.airfoil.chord;

wake = build_steady_wake( ...
    body, ...
    flow, ...
    wakeCfg);


%% ========================================================================
%  5) PRESCRIBE BODY SOURCE STRENGTHS
% =========================================================================

% Current stationary-body boundary condition:
%
%       sigma_i = -U_inf dot n_i
%
% This routine is intentionally isolated because it will later become:
%
%       sigma_i = -(U_inf - V_body,i) dot n_i
%
sigma = source_strengths( ...
    body, ...
    flow);


%% ========================================================================
%  6) BODY-BODY INFLUENCE COEFFICIENTS
% =========================================================================

% Influence of body source/doublet panels on body collocation points.
inflBB = body_influence(body);


%% ========================================================================
%  7) WAKE-BODY INFLUENCE COEFFICIENTS
% =========================================================================

% Influence of prescribed wake doublet panels on body collocation points.
inflBW = wake_influence( ...
    body, ...
    wake);


%% ========================================================================
%  8) ASSEMBLE DIRICHLET SYSTEM
% =========================================================================

% Assemble:
%
%       A * mu = b
%
% Unknowns:
%
%       mu(1:N) = body-panel doublet strengths
%       mu(end) = common steady-wake doublet strength
%
% Final row imposes the trailing-edge Kutta condition.
[A, b] = assemble_dirichlet_system( ...
    inflBB, ...
    inflBW, ...
    sigma);


%% ========================================================================
%  9) SOLVE FOR DOUBLET STRENGTHS
% =========================================================================

mu = A \ b;


%% ========================================================================
%  10) AERODYNAMIC POSTPROCESSING
% =========================================================================

aero = postprocess_aero( ...
    body, ...
    mu, ...
    sigma, ...
    inflBB, ...
    inflBW, ...
    flow, ...
    cfg);


%% ========================================================================
%  11) PACKAGE RESULTS
% =========================================================================

% ----- Original case definition -----

result.cfg = cfg;


% ----- Derived flow information -----

result.flow = flow;


% ----- Geometry -----

result.body = body;

result.wake = wake;


% ----- Singularity strengths -----

result.sigma = sigma;

result.mu = mu;

result.muBody = ...
    mu(1:body.nPanels);

result.muWake = ...
    mu(end);


% ----- Linear system -----

result.system.A = A;

result.system.b = b;


% ----- Aerodynamic data -----

result.aero = aero;


% Frequently used quantities promoted to top level
result.Cp = aero.Cp;

result.Vt = aero.Vt;

result.Gamma = aero.Gamma;

result.CL_gamma = aero.CL_gamma;

result.CL_pressure = aero.CL_pressure;

result.CD = aero.CD_pressure;

result.CM = aero.CM_pressure;


% ----- Cp plotting coordinates -----

% Cp is evaluated at body-panel collocation points.
result.xCp = ...
    body.cp(:,1) ./ cfg.airfoil.chord;


end