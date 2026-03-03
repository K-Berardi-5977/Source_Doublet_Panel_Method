% === INPUTS ===
% mu - vector of doublet strengths (including wake)
% sigma - body panel source strengths (numPan x 1)
% infl - struct containing influence coefficient data and local panel geometry
% mesh - struct containing panel geometry data 
% numPan - number of body panels
% num_d - number of doublets (body + wake)
% sys - struct containing matrix that entirely defines the zero internal perturbation potential Dirichlet condition matrix equations 

function aero = postprocess_aero(mu, sigma, infl, mesh, U, alphaR, numPan, num_d, sys)
% Compute wake induced surface velocity in receiver panel coordinate frame

ux_w = (mu(num_d)/(2*pi)).*(sys.zw./sys.rw2);
uz_w = -(mu(num_d)/(2*pi)).*(sys.xw./sys.rw2);
Vt_w = ux_w .* cos(mesh.th) + uz_w .* sin(mesh.th);

% Compute body doublet contribution to surface velocity
Vt_d = infl.BtJ * mu(1:num_d-1);

% Compute source contribution to surface velocity
Vt_s = infl.BtG * sigma;

% Compute free-stream contribution to surface velocity
U_s  = U .* cos(mesh.th - alphaR);

% Numerical differentiation correction term
mu_b = mu(1:numPan); % Isolate body doublets
g = zeros(numPan,1); % Vector to store jump terms
g(2:numPan-1) = (mu_b(3:end) - mu_b(1:end-2)) ./ (infl.Lp(2:end-1) + infl.Lp(1:end-2));
g(1) = (mu_b(2)-mu_b(1))/infl.Lp(1);
g(end) = (mu_b(end)-mu_b(end-1))/infl.Lp(end);

% Compute surface velocity
% Total tangential velocity
Vt = Vt_d + Vt_w + Vt_s + U_s + 0.5*g;

% Pressure coefficient
Cp = 1 - (Vt./U).^2;

% Circulation and lift coefficient
gamma = -sum(Vt .* infl.Lp(1:numPan));
cl    = 2*gamma;

% Package velocity, circulation, and aerodynamic load data
aero.Vt    = Vt;
aero.Cp    = Cp;
aero.gamma = gamma;
aero.cl    = cl;
end