function inflBW = wake_influence(body, wake)
%WAKE_INFLUENCE
% Compute influence of wake doublet panels on body collocation points.
%
% INPUTS
%   body : body panel geometry from panel_geometry()
%   wake : wake panel geometry from panel_geometry()
%
% OUTPUT
%   inflBW.A
%       Wake doublet potential influence at body collocation points
%
%   inflBW.BtJ
%       Wake-induced velocity resolved along body-panel tangent
%
%   inflBW.BnJ
%       Wake-induced velocity resolved along body-panel normal
%
% Unlike body_influence(), there are:
%   - no wake source terms
%   - no self-influence corrections
%
% because the wake panels and body receiver panels are distinct.


%% ========================================================================
%  1) BODY RECEIVER POINTS RELATIVE TO WAKE PANEL START
% =========================================================================

% Rows    -> body receiver collocation point i
% Columns -> wake observer panel j

Xg = body.cp(:,1) - wake.p1(:,1).';
Zg = body.cp(:,2) - wake.p1(:,2).';


%% ========================================================================
%  2) ROTATE INTO WAKE-PANEL LOCAL COORDINATES
% =========================================================================

Xp = ...
      Xg .* wake.cosTheta.' ...
    + Zg .* wake.sinTheta.';

Zp = ...
     -Xg .* wake.sinTheta.' ...
    + Zg .* wake.cosTheta.';


%% ========================================================================
%  3) WAKE PANEL LENGTH
% =========================================================================

X2p = wake.length.';


%% ========================================================================
%  4) DISTANCES TO WAKE PANEL ENDPOINTS
% =========================================================================

r1 = sqrt( ...
    Xp.^2 + Zp.^2);

r2 = sqrt( ...
    (Xp - X2p).^2 + Zp.^2);


%% ========================================================================
%  5) ANGLES TO WAKE PANEL ENDPOINTS
% =========================================================================

theta1 = atan2(Zp, Xp);

theta2 = atan2(Zp, Xp - X2p);


%% ========================================================================
%  6) WAKE DOUBLET POTENTIAL INFLUENCE
% =========================================================================

A = -(1/(2*pi)) .* ...
    (theta2 - theta1);


%% ========================================================================
%  7) WAKE DOUBLET VELOCITY IN WAKE LOCAL FRAME
% =========================================================================

J = (1/(2*pi)) .* ( ...
      Zp ./ r1.^2 ...
    - Zp ./ r2.^2 );


K = -(1/(2*pi)) .* ( ...
      Xp ./ r1.^2 ...
    - (Xp - X2p) ./ r2.^2 );


%% ========================================================================
%  8) ROTATE WAKE-INDUCED VELOCITY INTO BODY PANEL FRAME
% =========================================================================

% Difference between body receiver-panel angle and wake observer-panel angle
dTheta = ...
    body.theta - wake.theta.';

cij = cos(dTheta);
sij = sin(dTheta);


% Tangential velocity along body panel
BtJ = ...
    J .* cij ...
    + K .* sij;


% Normal velocity relative to body panel
BnJ = ...
   -J .* sij ...
    + K .* cij;


%% ========================================================================
%  9) PACKAGE OUTPUTS
% =========================================================================

inflBW.A = A;

inflBW.BtJ = BtJ;

inflBW.BnJ = BnJ;


end