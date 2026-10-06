function wake = build_steady_wake(body, flow, wakeCfg)
%BUILD_STEADY_WAKE Construct the prescribed straight steady wake.
%
% INPUTS
%   body
%       Body panel-geometry struct produced by panel_geometry()
%
%   flow
%       Flow struct containing:
%           flow.alphaRad
%
%   wakeCfg
%       Wake configuration containing:
%           wakeCfg.nPanels
%           wakeCfg.length
%
%
% OUTPUT
%   wake
%       Wake panel-geometry struct produced by panel_geometry()
%
%
% FUNCTION ROLE
%
%   This function determines WHERE the wake boundary points are.
%
%   It does not independently recompute panel lengths, angles,
%   collocation points, tangents, etc.
%
%   Those generic geometry operations are handled by:
%
%       panel_geometry()


%% ========================================================================
%  1) WAKE DIRECTION
% =========================================================================

% Current steady formulation assumes the wake extends directly along
% the freestream direction.

wakeDirection = [ ...
    cos(flow.alphaRad), ...
    sin(flow.alphaRad)];


%% ========================================================================
%  2) WAKE STARTING POINT
% =========================================================================

% Wake begins at the trailing edge of the body.
%
% With the current clockwise airfoil ordering, the final body boundary
% point is the trailing edge.

TE = body.Bp(end,:);


%% ========================================================================
%  3) WAKE BOUNDARY-POINT LOCATIONS
% =========================================================================

% Distance from the trailing edge to each wake boundary point.
%
% Produces Nw+1 points for Nw wake panels.

s = linspace( ...
    0, ...
    wakeCfg.length, ...
    wakeCfg.nPanels + 1).';


% Global wake boundary-point coordinates.
%
% Implicit expansion:
%
%       (Nw+1 x 1) .* (1 x 2)
%
% gives:
%
%       (Nw+1 x 2)

BpWake = TE + s .* wakeDirection;


%% ========================================================================
%  4) CONVERT WAKE POINTS INTO PANEL GEOMETRY
% =========================================================================

% panel_geometry() now performs all generic panel calculations:
%
%   wake.Bp
%   wake.p1
%   wake.p2
%   wake.cp
%   wake.length
%   wake.theta
%   wake.cosTheta
%   wake.sinTheta
%   wake.tangent
%   wake.normal
%   wake.nPanels

wake = panel_geometry(BpWake);


end