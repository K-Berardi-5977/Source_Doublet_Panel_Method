function inflBW = SDPM_Influence_BodyWake(mesh, wake)
% Influence of wake doublet panels (wake) on body collocation points/panels (mesh).
%
% Outputs (all size numPan x Nw):
%   inflBW.A    : doublet potential influence
%   inflBW.BtJ  : induced tangential velocity coefficient in body-panel tangential direction
%   inflBW.BnJ  : induced normal velocity coefficient in body-panel normal direction (optional)
%
% Assumes both mesh and wake were created by build_panels(), so they contain:
%   .th, .cos_th, .sin_th, .Bp1_x, .Bp1_z, .Bp2_x, .Bp2_z, .cx, .cz, .co

    numPan = size(mesh.co, 1);   % number of body panels (receivers)
    Nw     = size(wake.co, 1);   % number of wake panels (observers)

    % --- Transform body collocation points to wake-panel local coordinates ---
    % Offsets from wake panel first endpoint to body collocation point (global frame)
    Xg = mesh.cx - wake.Bp1_x.';     % (numPan x Nw)
    Zg = mesh.cz - wake.Bp1_z.';     % (numPan x Nw)

    % Wake panel vector components (global)
    X2g = wake.Bp2_x - wake.Bp1_x;   % (Nw x 1)
    Z2g = wake.Bp2_z - wake.Bp1_z;   % (Nw x 1)

    % Wake panel length projected onto wake local x-axis (should be positive)
    X2p_vec = X2g .* wake.cos_th + Z2g .* wake.sin_th;   % (Nw x 1)

    % Rotate global offsets into wake local frame
    Xp = Xg .* wake.cos_th.' + Zg .* wake.sin_th.';      % (numPan x Nw)
    Zp = -Xg .* wake.sin_th.' + Zg .* wake.cos_th.';     % (numPan x Nw)

    % Expand wake panel length across body rows
    X2p = ones(numPan,1) * X2p_vec.';                    % (numPan x Nw)

    % --- Constant-strength doublet integrals in wake-local coords ---
    rw1  = sqrt(Xp.^2 + Zp.^2);
    rw2  = sqrt((Xp - X2p).^2 + Zp.^2);

    thpw  = atan2(Zp, Xp);
    thpw2 = atan2(Zp, Xp - X2p);

    % Doublet potential influence
    A = -(1/(2*pi)) * (thpw2 - thpw);                      % (numPan x Nw)

    % Velocity influence in wake-local coordinates
    J = (1/(2*pi)) .* ( (Zp./(rw1.^2)) - (Zp./(rw2.^2)) ); % tangential (wake-local)
    K = -(1/(2*pi)) .* ( (Xp./(rw1.^2)) - ((Xp - X2p)./(rw2.^2)) ); % normal (wake-local)

    % --- Rotate induced velocity into body-panel coordinates ---
    % Need cos/sin of (theta_body - theta_wake)
    dth = mesh.th - wake.th_w.';                           % (numPan x Nw)
    c = cos(dth);
    s = sin(dth);

    BtJ = J .* c + K .* s;                               % tangential along body panel
    BnJ = -J .* s + K .* c;                              % normal along body panel

    inflBW.A   = A;
    inflBW.BtJ = BtJ;
    inflBW.BnJ = BnJ;
end