% ===== SDPM AoA SWEEP: Cp OVERLAY ===== %

clc
clear

%% ===== BASE CONFIGURATION =====

cfg = default_sdpm_config();

% Modify any common case settings here
cfg.airfoil.code      = '0012';
cfg.airfoil.chord     = 1;
cfg.airfoil.nPanels   = 90;
cfg.airfoil.spacing   = 'cosine';
cfg.airfoil.closedTE  = true;

cfg.flow.U            = 1;
cfg.flow.rho          = 1.225;

cfg.wake.nPanels      = 20;
cfg.wake.lengthChord  = 5;


%% ===== ANGLE-OF-ATTACK VECTOR =====

alphaVec = [-5 0 5 10];


%% ===== RUN SWEEP =====

results = cell(size(alphaVec));

figure(1)
clf
hold on

set(gca,'YDir','reverse')

for k = 1:numel(alphaVec)

    % Update only angle of attack
    cfg.flow.alphaDeg = alphaVec(k);

    % Run solver
    results{k} = solve_sdpm_steady(cfg);

    % Overlay Cp distribution
    plot( ...
        results{k}.xCp, ...
        results{k}.Cp, ...
        'LineWidth', 1.2, ...
        'DisplayName', ...
        sprintf('\\alpha = %.1f^\\circ', alphaVec(k)));

end


%% ===== PLOT FORMATTING =====

xlabel('x/c')
ylabel('C_p')

title(['NACA ', cfg.airfoil.code, ...
    ' Pressure Distributions'])

legend('Location','best')

grid on
axis padded

hold off