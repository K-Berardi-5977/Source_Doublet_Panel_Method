function cfg = prompt_sdpm_config()
%PROMPT_SDPM_CONFIG Interactively configure an SDPM simulation case.
%
% OUTPUT
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
% The function begins with default_sdpm_config(), then allows the user
% to override individual values. Pressing ENTER retains the default.


%% ========================================================================
%  LOAD DEFAULT CONFIGURATION
% =========================================================================

cfg = default_sdpm_config();

fprintf('\n');
fprintf('====================================================\n');
fprintf('             SDPM CASE CONFIGURATION\n');
fprintf('====================================================\n');
fprintf('Press ENTER to retain the value shown in brackets.\n');
fprintf('====================================================\n\n');


%% ========================================================================
%  AIRFOIL GEOMETRY
% =========================================================================

fprintf('--- AIRFOIL GEOMETRY ---\n\n');


% ----- NACA 4-digit designation -----

userInput = input( ...
    sprintf('NACA 4-digit code [%s]: ', cfg.airfoil.code), ...
    's');

if ~isempty(userInput)

    userInput = strtrim(userInput);

    if length(userInput) ~= 4 || ...
            any(~isstrprop(userInput, 'digit'))

        error('NACA designation must contain exactly four digits.');

    end

    cfg.airfoil.code = userInput;

end


% ----- Chord length -----

userInput = input( ...
    sprintf('Chord length [%.4g m]: ', cfg.airfoil.chord), ...
    's');

if ~isempty(userInput)

    value = str2double(userInput);

    if ~isfinite(value) || value <= 0
        error('Chord length must be a positive number.');
    end

    cfg.airfoil.chord = value;

end


% ----- Number of body panels -----

userInput = input( ...
    sprintf('Number of body panels [%d]: ', cfg.airfoil.nPanels), ...
    's');

if ~isempty(userInput)

    value = str2double(userInput);

    if ~isfinite(value) || ...
            value <= 0 || ...
            mod(value,1) ~= 0

        error('Number of body panels must be a positive integer.');

    end

    if mod(value,2) ~= 0
        error('Number of body panels must be even.');
    end

    cfg.airfoil.nPanels = value;

end


% ----- Panel spacing -----

userInput = input( ...
    sprintf('Panel spacing: cosine or uniform [%s]: ', ...
    cfg.airfoil.spacing), ...
    's');

if ~isempty(userInput)

    userInput = lower(strtrim(userInput));

    if ~ismember(userInput, {'cosine','uniform'})
        error('Panel spacing must be ''cosine'' or ''uniform''.');
    end

    cfg.airfoil.spacing = userInput;

end


% ----- Closed trailing edge -----

defaultTE = double(cfg.airfoil.closedTE);

userInput = input( ...
    sprintf('Closed trailing edge? 1 = yes, 0 = no [%d]: ', ...
    defaultTE), ...
    's');

if ~isempty(userInput)

    value = str2double(userInput);

    if ~ismember(value, [0 1])
        error('Trailing-edge selection must be either 0 or 1.');
    end

    cfg.airfoil.closedTE = logical(value);

end


%% ========================================================================
%  FLOW CONDITIONS
% =========================================================================

fprintf('\n--- FLOW CONDITIONS ---\n\n');


% ----- Angle of attack -----

userInput = input( ...
    sprintf('Angle of attack [%.4g deg]: ', cfg.flow.alphaDeg), ...
    's');

if ~isempty(userInput)

    value = str2double(userInput);

    if ~isfinite(value)
        error('Angle of attack must be numeric.');
    end

    cfg.flow.alphaDeg = value;

end


% ----- Freestream velocity -----

userInput = input( ...
    sprintf('Freestream velocity [%.4g m/s]: ', cfg.flow.U), ...
    's');

if ~isempty(userInput)

    value = str2double(userInput);

    if ~isfinite(value) || value <= 0
        error('Freestream velocity must be positive.');
    end

    cfg.flow.U = value;

end


% ----- Fluid density -----

userInput = input( ...
    sprintf('Fluid density [%.4g kg/m^3]: ', cfg.flow.rho), ...
    's');

if ~isempty(userInput)

    value = str2double(userInput);

    if ~isfinite(value) || value <= 0
        error('Fluid density must be positive.');
    end

    cfg.flow.rho = value;

end


%% ========================================================================
%  WAKE CONFIGURATION
% =========================================================================

fprintf('\n--- WAKE CONFIGURATION ---\n\n');


% ----- Number of wake panels -----

userInput = input( ...
    sprintf('Number of wake panels [%d]: ', cfg.wake.nPanels), ...
    's');

if ~isempty(userInput)

    value = str2double(userInput);

    if ~isfinite(value) || ...
            value < 1 || ...
            mod(value,1) ~= 0

        error('Number of wake panels must be a positive integer.');

    end

    cfg.wake.nPanels = value;

end


% ----- Wake length -----

userInput = input( ...
    sprintf('Wake length in chord lengths [%.4g c]: ', ...
    cfg.wake.lengthChord), ...
    's');

if ~isempty(userInput)

    value = str2double(userInput);

    if ~isfinite(value) || value <= 0
        error('Wake length must be positive.');
    end

    cfg.wake.lengthChord = value;

end


%% ========================================================================
%  AERODYNAMIC REFERENCE VALUES
% =========================================================================

fprintf('\n--- AERODYNAMIC REFERENCE VALUES ---\n\n');


% ----- Moment reference x-location -----

userInput = input( ...
    sprintf('Moment reference x-location [%.4g c]: ', ...
    cfg.reference.xMomentChord), ...
    's');

if ~isempty(userInput)

    value = str2double(userInput);

    if ~isfinite(value)
        error('Moment reference x-location must be numeric.');
    end

    cfg.reference.xMomentChord = value;

end


% ----- Moment reference z-location -----

userInput = input( ...
    sprintf('Moment reference z-location [%.4g c]: ', ...
    cfg.reference.zMomentChord), ...
    's');

if ~isempty(userInput)

    value = str2double(userInput);

    if ~isfinite(value)
        error('Moment reference z-location must be numeric.');
    end

    cfg.reference.zMomentChord = value;

end


%% ========================================================================
%  DISPLAY FINAL CONFIGURATION
% =========================================================================

fprintf('\n');
fprintf('====================================================\n');
fprintf('              FINAL SDPM CONFIGURATION\n');
fprintf('====================================================\n');

fprintf('Airfoil:              NACA %s\n', ...
    cfg.airfoil.code);

fprintf('Chord:                %.4g m\n', ...
    cfg.airfoil.chord);

fprintf('Body panels:          %d\n', ...
    cfg.airfoil.nPanels);

fprintf('Panel spacing:        %s\n', ...
    cfg.airfoil.spacing);

fprintf('Closed trailing edge: %d\n', ...
    cfg.airfoil.closedTE);

fprintf('\n');

fprintf('Angle of attack:      %.4g deg\n', ...
    cfg.flow.alphaDeg);

fprintf('Freestream velocity:  %.4g m/s\n', ...
    cfg.flow.U);

fprintf('Fluid density:        %.4g kg/m^3\n', ...
    cfg.flow.rho);

fprintf('\n');

fprintf('Wake panels:          %d\n', ...
    cfg.wake.nPanels);

fprintf('Wake length:          %.4g c\n', ...
    cfg.wake.lengthChord);

fprintf('\n');

fprintf('Moment x-reference:   %.4g c\n', ...
    cfg.reference.xMomentChord);

fprintf('Moment z-reference:   %.4g c\n', ...
    cfg.reference.zMomentChord);

fprintf('====================================================\n\n');

end