% === INPUTS ===
% mu - vector of doublet strengths (including wake)
% sigma - body panel source strengths (numPan x 1)
% infl - struct containing influence coefficient data and local panel geometry
% mesh - struct containing panel geometry data 
% numPan - number of body panels
% num_d - number of doublets (body + wake)
% sys - struct containing matrix that entirely defines the zero internal perturbation potential Dirichlet condition matrix equations 

function aero = postprocess_aero(mu, sigma, inflBB, inflBW, mesh, U, alphaR, numPan, num_d, sys, rho)
%% ===== Circulation Based Aerodynamic Loads =====

% Establish quarter-chord reference point for computation of moment
x_ref = 0.25;  
z_ref = 0.0;

% Compute wake induced surface velocity in receiver panel coordinate frame
Vt_w = sum(inflBW.BtJ, 2)*mu(num_d)

% Compute body doublet contribution to surface velocity
Vt_d = inflBB.BtJ * mu(1:num_d-1);

% Compute source contribution to surface velocity
Vt_s = inflBB.BtG * sigma;

% Compute free-stream contribution to surface velocity
U_s  = U .* cos(mesh.th - alphaR);

% Numerical differentiation correction term
mu_b = mu(1:numPan); % Isolate body doublets
dmu_dx = zeros(numPan,1); % Vector to store jump terms
dmu_dx(2:numPan-1) = (mu_b(3:end) - mu_b(1:end-2)) ./ (inflBB.Lp(2:end-1) + inflBB.Lp(1:end-2));
dmu_dx(1) = (mu_b(2)-mu_b(1))/inflBB.Lp(1);
dmu_dx(end) = (mu_b(end)-mu_b(end-1))/inflBB.Lp(end);

% Compute surface velocity
% Total tangential velocity
Vt = Vt_d + Vt_w + Vt_s + U_s + 0.5*dmu_dx;

% Pressure coefficient
Cp = 1 - (Vt./U).^2;

% Circulation and lift coefficient
gamma = sum(Vt .* inflBB.Lp(1:numPan));
cl    = 2*gamma;

% Package velocity, circulation, and aerodynamic non-dimensional coefficient data
aero.Vt    = Vt;
aero.Cp    = Cp;
aero.gamma = gamma;
aero.cl    = cl;

%% ===== Pressure force integration aerodynamic loads =====

% Dynamic pressure of free-stream
q = 0.5 * rho * U^2;

% Panel lengths
ds = inflBB.Lp(1:numPan);

% Outward unit normal for clockwise panel ordering
nx =  -sin(mesh.th);
nz = cos(mesh.th);
n = [nx, nz]; % Unit normal vector for each panel

% Split Cp into constant and velocity-dependent parts for debugging
Cp_const = ones(numPan,1);
Cp_dyn   = -(Vt./U).^2;

% Force contributions from each part
dFx_const = -q .* Cp_const .* nx .* ds;
dFz_const = -q .* Cp_const .* nz .* ds;

dFx_dyn   = -q .* Cp_dyn .* nx .* ds;
dFz_dyn   = -q .* Cp_dyn .* nz .* ds;

% Total panel forces
dFx = dFx_const + dFx_dyn;
dFz = dFz_const + dFz_dyn;
% Total force components in global x-z coordinates
Fx = sum(dFx);
Fz = sum(dFz);

% Resolve into drag and lift relative to freestream direction
D =  Fx*cos(alphaR) + Fz*sin(alphaR);
L = -Fx*sin(alphaR) + Fz*cos(alphaR);

Fx_const = sum(dFx_const);
Fz_const = sum(dFz_const);

% ----------------------------------------
% Moment about reference point (per unit span)
% ----------------------------------------
rx = mesh.cx - x_ref;
rz = mesh.cz - z_ref;

dM = rx .* dFz - rz .* dFx;
M  = sum(dM);

% ----------------------------------------
% Force coefficients (2-D: based on chord per unit span)
% ----------------------------------------
c_ref = max(mesh.Bp(:,1)) - min(mesh.Bp(:,1));   % simple chord estimate

cd = D / (q * c_ref);
cl_int = L / (q * c_ref)
cm = M / (q * c_ref^2);
cl_gamma = 2*gamma(U*c_ref)

% Package force and moment data
aero.dFx = dFx; % global x-direction force-per-unit-span component at the collocation point of each panel
aero.dFz = dFz; % global z-direction force-per-unit-span component at the collocation point of each panel
aero.Fx = Fx; % net global x-direction force-per-unit-span at the centroid of the airfoil
aero.Fz = Fz; % net global z-direction force-per-unit-span at the centroid of the airfoil
aero.D = D; % net drag-per-unit-span (should be >= 1e^-3 for inviscid)
aero.L = L; % net lift-per-unit-span 
aero.M = M; % net quater chord airfoil moment 
aero.Fx_const = Fx_const;
aero.Fz_const = Fz_const;
aero.cm = cm; % Quarter chord moment coefficient

end