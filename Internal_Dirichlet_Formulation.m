% === INPUTS ===
% infl - struct containing influence coefficient data and local panel geometry
% sigma - body panel source strengths (numPan x 1)
% mesh - struct containing panel geometry data 
% numPan - number of body panels
% num_d - number of doublets (body + wake)

% === OUTPUTS ===
% sys - struct containing matrix that entirely defines the zero internal perturbation potential Dirichlet condition matrix equations for potential, and contains wake panel coordinates 
% sys.A_full - matrix containing influence coefficients (body sources/doublets, and wake) with Kutta condition applied
% sys.xw  - wake x-distance from body collocation points (numPan x 1)
% sys.zw  - wake z-distance from body collocation points (numPan x 1)
% sys.rw - wake distance from body collocation points


function sys = Internal_Dirichlet_Formulation(infl, sigma, mesh, numPan, num_d)

% Construct right-hand-side source terms
RHS = infl.SourceInfluence * sigma;

% Compute distance between wake panel collocation points and body panel
% collocation points
xw = mesh.cx - mesh.Bp2(numPan,1);
zw = mesh.cz - mesh.Bp2(numPan,2);
rw2 = xw.^2 + zw.^2;
thw = -atan(zw./xw); % Angle between wake and receiver panel

% Create A_full to store matrix eqn data
A_full = zeros(numPan+1, numPan+1);

A_full(1:numPan,1:numPan) = infl.A; % Doublet body panel influence coefficients
A_full(1:numPan,num_d) = -(1/(2*pi))*thw; % Wake panel influence coefficients
A_full(1:numPan,num_d+1) = RHS; % Source influence (strength x influence coefficient)

% Apply Kutta Condition to trailing edge
A_full(numPan+1,:) = 0;
A_full(numPan+1,1) = -1;
A_full(numPan+1,numPan) = 1;
A_full(numPan+1,num_d) = -1;

% Package A_full and necessary wake geometry definitions

sys.A_full = A_full;
sys.xw  = xw;
sys.zw  = zw;
sys.rw2 = rw2;

end