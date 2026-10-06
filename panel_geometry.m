function panels = panel_geometry(Bp)
%PANEL_GEOMETRY Construct reusable geometry for a connected panel set.
%
% INPUT
%   Bp : (N+1 x 2) array of panel boundary points
%
%        Bp(:,1) = x-coordinate
%        Bp(:,2) = z-coordinate
%
% OUTPUT
%   panels : struct containing panel geometry
%
%       panels.Bp
%           Boundary points
%
%       panels.p1
%           First endpoint of each panel
%
%       panels.p2
%           Second endpoint of each panel
%
%       panels.cp
%           Panel collocation points
%
%       panels.length
%           Panel lengths
%
%       panels.theta
%           Panel angle measured from +x axis
%
%       panels.cosTheta
%       panels.sinTheta
%           Reusable trigonometric quantities
%
%       panels.tangent
%           Unit tangent vectors
%
%       panels.normal
%           Unit normal vectors
%
%       panels.nPanels
%           Number of panels
%
%
% IMPORTANT
%   The normal convention used here is
%
%       n = [-sin(theta), cos(theta)]
%
%   which is the outward normal for the CLOCKWISE body-panel ordering
%   currently used by the SDPM airfoil geometry.


%% ========================================================================
%  BASIC VALIDATION
% =========================================================================

if size(Bp,2) ~= 2
    error('Bp must be an (N+1) x 2 array of [x,z] coordinates.');
end

if size(Bp,1) < 2
    error('At least two boundary points are required.');
end


%% ========================================================================
%  PANEL ENDPOINTS
% =========================================================================

panels.Bp = Bp;

panels.p1 = Bp(1:end-1,:);

panels.p2 = Bp(2:end,:);


%% ========================================================================
%  PANEL VECTORS
% =========================================================================

dP = panels.p2 - panels.p1;

dx = dP(:,1);

dz = dP(:,2);


%% ========================================================================
%  PANEL LENGTHS
% =========================================================================

panels.length = hypot(dx,dz);


% Prevent degenerate panels
if any(panels.length <= eps)
    error('Panel geometry contains one or more zero-length panels.');
end


%% ========================================================================
%  PANEL ANGLES
% =========================================================================

panels.theta = atan2(dz,dx);

panels.cosTheta = cos(panels.theta);

panels.sinTheta = sin(panels.theta);


%% ========================================================================
%  UNIT TANGENT VECTORS
% =========================================================================

panels.tangent = [ ...
    panels.cosTheta, ...
    panels.sinTheta];


%% ========================================================================
%  UNIT NORMAL VECTORS
% =========================================================================

% For clockwise body-panel ordering:
%
%       tangent = [ cos(theta),  sin(theta)]
%
%       normal  = [-sin(theta),  cos(theta)]
%
panels.normal = [ ...
    -panels.sinTheta, ...
     panels.cosTheta];


%% ========================================================================
%  COLLOCATION POINTS
% =========================================================================

panels.cp = ...
    0.5 .* (panels.p1 + panels.p2);


%% ========================================================================
%  PANEL COUNT
% =========================================================================

panels.nPanels = ...
    size(panels.p1,1);


end