function [Nx , Ny, Vxy, theta_pj, psi, THETA, Cpxy_mask] = PM_streamlines(xi, yi, Xj, Yj, phiR, S, lambda, U, alpha, Cp, numPan)
alphaR = alpha*(pi/180); %convert angle of attack to rads
[Nx , Ny] = SPM_xyCoeff(xi, yi, Xj, Yj, phiR, S); %Function to Compute X and Y velocity influence coefficients for Meshgrid

nGridX = numPan;                                                        % X-grid for streamlines and contours
nGridY = numPan;                                                        % Y-grid for streamlines and contours
xVals  = [-0.3; 1.3];                                                   % X-grid extents [min, max]
yVals  = [-.8; 0.8];                                                    % Y-grid extents [min, max]

% Generate the grid points
Xgrid   = linspace(xVals(1),xVals(2),nGridX)';                          % X-values in evenly spaced grid
Ygrid   = linspace(yVals(1),yVals(2),nGridY)';                          % Y-values in evenly spaced grid
[XX,YY] = meshgrid(Xgrid,Ygrid);                                        % Create meshgrid from X and Y grid arrays
rr = sqrt(XX.^2 + YY.^2);
THETA = atan2d(YY, XX);
% Streamline parameters
stepsize = 0.01;                                                        % Step size for streamline propagation
maxVert  = nGridX*nGridY*10;                                            % Maximum vertices
slPct    = 30;                                                          % Percentage of streamlines of the grid
Ysl      = linspace(yVals(1),yVals(2),floor((slPct/100)*nGridY))';      % Create array of Y streamline starting points



% Initialize velocities
Vx = zeros(nGridX,nGridY); % Initialize gridpoint x-velocity matrix
Vy = zeros(nGridX,nGridY); % Initialize gridpoint y-velocity matrix

for m = 1:1:nGridX   %iterating over the jth control point
    for n = 1:1:nGridY  %for each control point, iterate over j=1:n panels
            XP = XX(m,n);   %Current iteration's x grid point
            YP = YY(m,n);   %Current iteration's y grid point
            [Nxx, Nyy] = SPM_xyCoeff(XP, YP, Xj, Yj, phiR,S); %calculate cartesian influence coefficients for each source at the current gridpoint
            [in,on] = inpolygon(XP,YP,Xj, Yj);% See if points are in or on the airfoil
             if (in == 1 || on == 1)
                 Vx(m,n) = 0; % Set X-velocity equal to zero                                              
                 Vy(m,n) = 0; % Set Y-velocity equal to zero                                                
            else                                                            % If the grid point is outside the airfoil
                Vx(m,n) = U*cos(alpha) + sum((lambda*Nxx)/(2*pi));         % Compute X-velocity
                Vy(m,n) = U*sin(alpha) + sum((lambda*Nyy)/(2*pi));         % Compute Y-velocity
            end
          Vxy(m,n) =  sqrt(Vx(m,n)^2 + Vy(m,n)^2); %velocity magnitude mesh
          Cpxy(m,n) = 1-(Vxy(m,n)/U)^2; %pressure coefficient mesh
    end

end


Cpxy_mask = Cpxy;
[inC, onC] = inpolygon(XX, YY, Xj, Yj);
Cpxy_mask(inC) = NaN;



%========= Stream Function Expression ==========%

lambda_dS = lambda(:).*S(:); %determine strength of each vortex 
psi = U*(YY*cos(alphaR)-XX*sin(alphaR)); %initialize stream function variable and account for sine term

%Calculate the sum of the vortex contributions to the velocity potential

for j=1:numPan
    dxp = (XX-xi(j)); %x-distance between grid points and jth panel
    dyp = (YY-yi(j)); %y-distance between grid points and jth panel
    theta_pj = atan2(dyp, dxp); %angle between jth source and all the grid points with respect to +x-axis
    psi = psi + (lambda_dS(j)/(2*pi))*theta_pj; %add contribution of jth panel to stream function
end

[in on] = inpolygon(XX, YY, Xj, Yj); %find mesh values within foil body
psi_mask = psi;
psi_mask(in) = NaN; %disregard streamlines solved at mesh values that fall within the surface 


%========== Plot Streamlines and Stagnation Points =========%
[row, col] = find(Cp>0.99); %Find coordinates in field that coincide with Cp ~ 1
aoa_D = alpha %convert angle of attack to degrees for plotting

figure; hold on; box on
contour(XX, YY, psi_mask, 30, 'LineWidth', 0.9); %streamlines
fill(Xj, Yj, [0.15 0.15 0.15], 'EdgeColor', 'k'); % foil body
for i = 1:length(row)
    plot(Xj(row(i)), Yj(row(i)), 'ro', 'MarkerSize', 3, MarkerFaceColor='r')
end
axis equal tight
xlabel('x/c'); ylabel('y/c');
title(['VPM Streamlines ($\alpha = ', num2str(aoa_D), '$ deg)' ], 'Interpreter','latex')
legend("Streamlines",'Airfoil Body','Stagnation Points')
end

