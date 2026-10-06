function cfg = default_sdpm_config()
% =================================================================== %
%% ========== Airfoil Characteristics ========== %
% =================================================================== %
cfg.airfoil.code      = '0012'; % 4-Digit NACA Airfoil Specification -- NACA 0012 is a canonical test-case with well documented solutions and test data
cfg.airfoil.chord     = 1; % Chord length [units: meters]
cfg.airfoil.nPanels   = 90; % Number of body panels
cfg.airfoil.spacing   = 'cosine'; % Cosine spacing selected as default due to more stable handling of potential gradient near trailing edge
cfg.airfoil.closedTE  = true; % Default configuration will have an assumed closed trailing edge

% =================================================================== %
%% ========== Flow Characteristics ========== %
% =================================================================== %
cfg.flow.U            = 1; % Free-stream velocity 
cfg.flow.alphaDeg     = 5;
cfg.flow.rho          = 1.225;

% =================================================================== %
%% ========== Wake Characteristics ========== %
% =================================================================== %
cfg.wake.nPanels      = 20;
cfg.wake.lengthChord  = 5;

% =================================================================== %
%% ========== Reference Geometry ========== %
% =================================================================== %
cfg.reference.xMomentChord = 0.25;
cfg.reference.zMomentChord = 0;

end