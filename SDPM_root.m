% ===== SOURCE-DOUBLET PANEL METHOD EXECUTION / VALIDATION SCRIPT ===== %
% This script defines a CASE, calls the steady SDPM solver, and performs
% case-specific plotting / comparison against experimental data.
%
% CHANGED:
%   - All solver inputs now live in cfg.
%   - Bypass interactive prompt with default configuratin
%   - Plotting remains here rather than inside geometry/solver functions.
clc; clear; 

choice = questdlg( ...
    'How would you like to Configure the SDPM Solver?', ...
    'SDPM Setup', ...
    'Use Defaults', ...
    'Enter Inputs', ...
    'Cancel', ...
    'Use Defaults')

switch choice
    case 'Use Defaults'
        cfg = default_sdpm_config(); % Default simulation case (see comments inside function for details)
    case 'Enter Inputs'
        cfg = prompt_sdpm_config(); % User defines the airfoil geometry, initial wake, and initial flow field characteristics
    otherwise
        disp('SDPM execution cancelled')
        return
end

result = solve_sdpm_steady(cfg); 

%% ========================================================================
%  VALIDATION / SMOKE-TEST PLOTS
% =========================================================================

%% ----- 1) AIRFOIL / PANEL GEOMETRY -----

figure(1)
clf
hold on

plot( ...
    result.body.Bp(:,1) ./ cfg.airfoil.chord, ...
    result.body.Bp(:,2) ./ cfg.airfoil.chord, ...
    '-ok', ...
    'LineWidth', 1.0, ...
    'MarkerSize', 3);

plot( ...
    result.body.cp(:,1) ./ cfg.airfoil.chord, ...
    result.body.cp(:,2) ./ cfg.airfoil.chord, ...
    '.r', ...
    'MarkerSize', 10);

axis equal
grid on
axis padded

xlabel('x/c')
ylabel('z/c')

title(['NACA ', cfg.airfoil.code, ' Panel Geometry'])

legend( ...
    'Panel boundary', ...
    'Collocation points', ...
    'Location', 'best');

hold off


%% ----- 2) PRESSURE COEFFICIENT DISTRIBUTION -----

figure(2)
clf
hold on

% Conventional aerodynamic Cp plotting direction
set(gca, 'YDir', 'reverse')

plot( ...
    result.xCp, ...
    result.Cp, ...
    '-ob', ...
    'LineWidth', 1.1, ...
    'MarkerSize', 4);

grid on
axis padded

xlabel('x/c')
ylabel('C_p')

title( ...
    ['NACA ', cfg.airfoil.code, ...
     ', \alpha = ', num2str(cfg.flow.alphaDeg), '^\circ']);

hold off


%% ----- 3) PRINT AERODYNAMIC VALIDATION QUANTITIES -----

fprintf('\n');
fprintf('============================================\n');
fprintf('           SDPM VALIDATION RESULTS\n');
fprintf('============================================\n');

fprintf('NACA:              %s\n', cfg.airfoil.code);
fprintf('Angle of attack:   %.4f deg\n', cfg.flow.alphaDeg);
fprintf('Body panels:       %d\n', result.body.nPanels);
fprintf('Wake panels:       %d\n', result.wake.nPanels);

fprintf('\n');

fprintf('Gamma:              %.8f\n', result.aero.Gamma);
fprintf('CL - circulation:   %.8f\n', result.aero.CL_gamma);
fprintf('CL - pressure:      %.8f\n', result.aero.CL_pressure);
fprintf('CD - pressure:      %.8e\n', result.aero.CD_pressure);
fprintf('CM - pressure:      %.8f\n', result.aero.CM_pressure);

fprintf('\n');

CL_difference = ...
    result.aero.CL_pressure - result.aero.CL_gamma;

CL_percent_difference = ...
    100 .* abs(CL_difference) ./ ...
    max(abs(result.aero.CL_gamma), eps);

fprintf('Delta CL:           %.8e\n', CL_difference);
fprintf('CL difference:      %.4f %%\n', CL_percent_difference);

fprintf('\n');

fprintf('Uniform-pressure force check:\n');
fprintf('Fx_const:           %.8e\n', result.aero.Fx_const);
fprintf('Fz_const:           %.8e\n', result.aero.Fz_const);

fprintf('============================================\n\n');

%% ====================================================================
%  FIGURE 3: Compare Cp Against Experimental Validation Data
% =====================================================================
%% ===== VALIDATION DATA =====

% Experimental / reference data format:
%   Column 1 = x/c
%   Column 2 = Cp

validationData = load('Cp_Gregory_Oreilly.dat');

figure(3)
clf
hold on

% Conventional Cp plotting direction
set(gca, 'YDir', 'reverse')


% ----- SDPM RESULT -----

plot( ...
    result.xCp, ...
    result.Cp, ...
    '-b', ...
    'LineWidth', 1.2, ...
    'DisplayName', 'SDPM');


% ----- INPUT VALIDATION DATA -----

plot( ...
    validationData(:,1), ...
    validationData(:,2), ...
    'ok', ...
    'MarkerSize', 5, ...
    'LineWidth', 1.0, ...
    'DisplayName', 'Validation Data');


%% ===== FORMATTING =====

xlabel('x/c')
ylabel('C_p')

title( ...
    ['NACA ', cfg.airfoil.code, ...
    ': SDPM vs Validation Data, \alpha = ', ...
    num2str(cfg.flow.alphaDeg), '^\circ'])

legend('Location', 'best')

grid on
axis padded

hold off