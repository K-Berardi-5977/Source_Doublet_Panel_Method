function inflBB = body_influence(body)
%BODY_INFLUENCE
% Compute body-panel source and doublet influence coefficients at body
% collocation points.
%
% INPUT
%   body : panel geometry struct from panel_geometry()
%
% OUTPUT
%   inflBB.A
%       Body doublet potential influence coefficients
%
%   inflBB.BtJ
%       Body doublet tangential velocity influence coefficients
%
%   inflBB.SourceInfluence
%       Body source potential influence coefficients
%
%   inflBB.BtG
%       Body source tangential velocity influence coefficients
%
% Self-influence terms are explicitly assigned because the receiver
% collocation point lies on its own source panel.


%% ========================================================================
%  1) BASIC DIMENSIONS
% =========================================================================

n = body.nPanels;

% Linear indices of diagonal entries in an n x n matrix
diagIdx = 1:(n+1):n^2;


%% ========================================================================
%  2) RECEIVER POINTS RELATIVE TO OBSERVER PANEL START
% =========================================================================

% Rows    -> receiver body collocation point i
% Columns -> observer body panel j

Xg = body.cp(:,1) - body.p1(:,1).';
Zg = body.cp(:,2) - body.p1(:,2).';


%% ========================================================================
%  3) ROTATE INTO OBSERVER-PANEL LOCAL COORDINATES
% =========================================================================

Xp = ...
      Xg .* body.cosTheta.' ...
    + Zg .* body.sinTheta.';

Zp = ...
     -Xg .* body.sinTheta.' ...
    + Zg .* body.cosTheta.';


%% ========================================================================
%  4) OBSERVER PANEL LENGTH
% =========================================================================

% body.length is n x 1.
% Transpose allows implicit expansion across all receiver rows.

X2p = body.length.';


%% ========================================================================
%  5) DISTANCES TO OBSERVER PANEL ENDPOINTS
% =========================================================================

r1 = sqrt( ...
    Xp.^2 + Zp.^2);

r2 = sqrt( ...
    (Xp - X2p).^2 + Zp.^2);


%% ========================================================================
%  6) ANGLES TO OBSERVER PANEL ENDPOINTS
% =========================================================================

theta1 = atan2(Zp, Xp);

theta2 = atan2(Zp, Xp - X2p);


%% ========================================================================
%  7) DOUBLET POTENTIAL INFLUENCE
% =========================================================================

A = -(1/(2*pi)) .* ...
    (theta2 - theta1);

% Body-panel doublet self influence
A(diagIdx) = 0.5;


%% ========================================================================
%  8) DOUBLET VELOCITY IN OBSERVER LOCAL FRAME
% =========================================================================

% Tangential component in observer-panel frame
J = (1/(2*pi)) .* ( ...
      Zp ./ r1.^2 ...
    - Zp ./ r2.^2 );

% Doublet tangential self influence
J(diagIdx) = 0;


% Normal component in observer-panel frame
K = -(1/(2*pi)) .* ( ...
      Xp ./ r1.^2 ...
    - (Xp - X2p) ./ r2.^2 );


% Preserve existing self-influence treatment
XpDiag = Xp(diagIdx);
XpDiag = XpDiag(:);

LpDiag = body.length(:);

K(diagIdx) = -(1/(2*pi)) .* ( ...
      1 ./ XpDiag ...
    - 1 ./ (XpDiag - LpDiag) );


%% ========================================================================
%  9) ROTATE DOUBLET VELOCITY INTO RECEIVER TANGENT DIRECTION
% =========================================================================

% theta_i - theta_j
dTheta = body.theta - body.theta.';

cij = cos(dTheta);
sij = sin(dTheta);

BtJ = ...
    J .* cij ...
    + K .* sij;


%% ========================================================================
%  10) SOURCE POTENTIAL INFLUENCE
% =========================================================================

SourceInfluence = (1/(2*pi)) .* ( ...
      Xp .* log(r1) ...
    - (Xp - X2p) .* log(r2) ...
    + Zp .* (theta2 - theta1) );


% Preserve existing source self-influence expression
SourceInfluence(diagIdx) = ...
    (1/pi) .* ...
    Xp(diagIdx) .* log(r1(diagIdx));


%% ========================================================================
%  11) SOURCE VELOCITY IN OBSERVER LOCAL FRAME
% =========================================================================

G = (1/(2*pi)) .* ...
    log(r1 ./ r2);


% Source tangential self influence
G(diagIdx) = ...
    (1/(2*pi)) .* ...
    log( ...
        XpDiag ./ ...
        abs(XpDiag - LpDiag));


H = (1/(2*pi)) .* ...
    (theta2 - theta1);

% Source normal self influence
H(diagIdx) = -0.5;


%% ========================================================================
%  12) ROTATE SOURCE VELOCITY INTO RECEIVER TANGENT DIRECTION
% =========================================================================

BtG = ...
    G .* cij ...
    + H .* sij;


%% ========================================================================
%  13) PACKAGE OUTPUTS
% =========================================================================

inflBB.A = A;

inflBB.BtJ = BtJ;

inflBB.SourceInfluence = SourceInfluence;

inflBB.BtG = BtG;


end