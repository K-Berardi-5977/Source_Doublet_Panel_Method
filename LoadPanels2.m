function [Xb, Yb, Xc, Yc, betaR, rot, S, numPan, n, t] = LoadPanels2(alpha)

%% ===== Initialize Parameters ===== %
load('foilData.dat'); % Load airfoil grid (boundary) points from data file
N = length(foilData(:,1)); % Number of boundary points
numPan = N-1; % Number of panels

Xb = foilData(:,1); % Global x-coordinates of grid points
Yb = foilData(:,2); % Global y-coordinates of grid points
Yb = flip(Yb); % Flip y-vector to accomodate cw iteration from trailing egde

Xc = zeros(numPan,1); % Initialize control points x-coordinate vector
Yc = zeros(numPan,1); % Initialize control points y-coordinate vector

%% ===== Generate Panel Geometry ===== %
dX = diff(Xb)'; % Compute panel length x-component 
dY = diff(Yb)'; % Compute panel length y-component 
S = hypot(dX, dY); % Compute panel length
t = ([dX; dY]./S)'; % Compute panel tangent vector
n = [-t(:,2) t(:,1)]; % Compute panel normal vector
beta = atan2d(dY, dX); % Compute angle (deg) between panel and global (positive) x-axis
beta = mod(beta, 360); % Make all angles positive and between 0-360 degrees
betaR = beta.*(pi/180); % Convert beta to radians

rot = [cos(betaR) sin(betaR); 
    -sin(betaR) cos(betaR)];  % Define rotation matrix to convert vector quantities to local panel coordinates



% Generate Control Points in Global Coordinates
for i = 1:numPan
    Xc(i) = 0.5*(Xb(i+1)+Xb(i)); % X-coordinate of control point of ith panel
    Yc(i) = 0.5*(Yb(i+1)+Yb(i)); % Y-coordinate of control point of ith panel
end




end