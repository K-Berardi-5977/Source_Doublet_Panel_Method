% === INPUTS ===
% mesh - struct containing panel geometry data 
% numPan - number of body panels
% diag_idx - diagonal index variable for any numPan x numPan matrix

% === OUTPUTS ===
% infl - struct containing influence coefficient data and local panel geometry
% infl.A - Body doublet potential influence coefficients, NOT WAKE POTENTIAL
% infl.BtJ - Body doublet tangent velocity influence coefficients
% infl.SourceInfluence = SourceInfluence;
% infl.BtG - Body source tangent velocity influence coefficients
% infl.Xp - Distance from first boundary point of observer panel to the control point of receiver panel in observer panel coordinates
% infl.Zp   = Vertical distance from observer panel to control point of receiver panel in observer panel coordinates;
% infl.X2p  = X2p; % Distance between boundary points of observer panel in observer panel coordinates
% infl.Lp   = Lp; % Panel lengths

function infl = SDPM_Influence(mesh, numPan, diag_idx)

% Transform collocation points to local panel coordinates
Xg = mesh.cx - mesh.Bp1_x.'; % (numPan x numPan)
Zg = mesh.cz - mesh.Bp1_z.'; % (numPan x numPan)

X2g = mesh.Bp2_x - mesh.Bp1_x; % (numPan x 1)
Z2g = mesh.Bp2_z - mesh.Bp1_z; % (numPan x 1)

% Perform coordinate rotation from global to local panel coordinates
Xp = Xg .* repmat(mesh.cos_th.', numPan,1) + Zg .*repmat(mesh.sin_th.', numPan,1);
Zp = -Xg .*repmat(mesh.sin_th.', numPan,1) + Zg .* repmat(mesh.cos_th.', numPan,1);

% Compute panel lengths in local panel coordinates (validation: should all be positive for proper panel orientation)
X2p_vec = X2g .* mesh.cos_th + Z2g .* mesh.sin_th; % (numPan × 1)
X2p = ones(numPan,1)*X2p_vec.'; % (numPan x numPan)

% Store panel lengths for future computations
Lp = X2p_vec(:);

% Distances between boundary points and collocation points (each column is jth panel, each row
% is ith collocation point)
r1 = sqrt(Xp.^2 + Zp.^2);
r2 = sqrt((Xp - X2p).^2 + Zp.^2);

% Angles between boundary points and collocation points (each column is jth panel, each row
% is ith collocation point)
thp = atan2(Zp, Xp);
thp2 = atan2(Zp, Xp - X2p);

% ===== Compute doublet potential influence coefficients in local panel coordinates =====
A = -(1/(2*pi))*(thp2 - thp);

A(diag_idx) = 0.5; % self influence coefficient

% ===== Compute doublet velocity influence coefficients in local panel coordinates =====
J = (1/(2*pi)).*((Zp./(r1.^2))-(Zp./(r2.^2))); % Surface velocity influence coefficient
J(diag_idx) = 0; % Surface velocity self influence

K = -(1/(2*pi)).*((Xp./(r1.^2))-((Xp-X2p)./(r2.^2))); % Normal velocity influence coefficient
K(diag_idx) = -(1/(2*pi)) * ( (1./Xp(diag_idx)) - (1./(Xp(diag_idx)-X2p(diag_idx))) ); % Normal velocity self influence

% === Rotate doublet velocity influence coefficients into receiver panel tangent direction ===
BtJ = J.*mesh.c + K.*mesh.s;

% ===== Compute source potential influence coefficients in local panel coordinates =====
SourceInfluence = (1/(2*pi))*( ...
    Xp .* log(r1) ...
    - (Xp - X2p) .* log(r2) ...
    + Zp .* (thp2 - thp) ...
    );

SourceInfluence(diag_idx) = (1/pi)*(Xp(diag_idx).*log(r1(diag_idx))); % self influence coefficient


% ===== Compute source velocity influence coefficients in local panel coordinates =====
G = (1/(2*pi)).*(log((r1./r2))); % Surface velocity influence coefficient
G(diag_idx) = (1/(2*pi)).*log(Xp(diag_idx)./(abs(Xp(diag_idx)-X2p(diag_idx)))); % Surface velocity self influence

H = (1/(2*pi))*(thp2-thp); % Normal velocity influence coefficient
H(diag_idx) = -0.5; % Normal velocity self influence

% === Rotate source velocity influence coefficients into receiver panel tangent direction ===
BtG = G.*mesh.c + H.*mesh.s;

% ===== Package influence coefficients and any necessary geometry updates into struct =====

% Influence coefficient exports
infl.A               = A;
infl.BtJ             = BtJ;
infl.SourceInfluence = SourceInfluence;
infl.BtG             = BtG;

% Geometry exports
infl.Xp   = Xp;
infl.Zp   = Zp;
infl.X2p  = X2p;
infl.Lp   = Lp;
end

