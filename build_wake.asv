% Function to build prescribed wake panel geometry

% === INPUTS ===
% mesh - struct containing panel geometry data 
% U -  free-stream velocity magnitude
% alphaR - pitch angle [radians]
% c - chord length
% Nw - number of wake panels
% Lw - wake length

% === Outputs ===
% wake - struct similar to 'mesh' containing prescribed wake geometry data
% wake.Bp1, wake.Bp2 - wake boundary points
% wake.th - angle between wake and +x axis
% wake.cos_th - cos(th) for each wake panel (useful for panel frame rotations) (Nw x 1)
% wake.sin_th - sin(th) for each wake panel (useful for panel frame rotations) (Nw x 1)
% wake.co - collocation point coordinates as (x,y) pairs (Nw x 2)
% wake.cx - collocation point x-coordinates (Nw x 1)
% wake.cz - collocation point z-coordinates (Nw x 1)

function wake = build_wake(mesh, U, alphaR, c, Nw, Lw)
% Initialize variables
u_hat = [cos(alphaR), sin(alphaR)]; % Free-stream direction, assume wake travels along free-stream
sw = (0:Nw)'*(Lw/Nw); % Wake panel lengths calculated as Lengthof wake/Number of wake panels
% Establish wake panel boundary points beginning with trailing edge of foil
TE = mesh.Bp2(end, :); % TE x,y coordinate pair -- validation: at t=0, TE = 1,0

Bp_w = TE + sw.*u_hat; % (Nw+1 x 2) matrix of wake panel boundary points


% Panel vectors and panel angles (global frame)
dxy_w     = diff(Bp_w);                          % (Nw   x 2)
th_w      = atan2(dxy_w(:,2), dxy_w(:,1));         % Panel angle wrt +x axis (numPan x 1)
wake.th_w = th_w;                                % Store panel angle
wake.cos_th = cos(th_w);                       % (Nw x 1) useful for rotating coordinate frame
wake.sin_th = sin(th_w);                       % (Nw x 1) useful for rotating coordinate frame

% Panel endpoint vectors
Bp_w1 = Bp_w(1:Nw, :);
Bp_w2 = Bp_w(2:Nw+1, :);
wake.Bp = Bp_w;                                % Boundary point coordinates as (x,z) pairs (NN-1 x 2)
wake.Bp1 = Bp_w1;                              % Store Bp1
wake.Bp2 = Bp_w2;                              % Store Bp2
wake.Bp1_x = Bp_w1(:,1);  wake.Bp1_z = Bp_w1(:,2); % Bp1 x and z coordinates (Nw x 1)
wake.Bp2_x = Bp_w2(:,1);  wake.Bp2_z = Bp_w2(:,2); % Bp2 x and z coordinates (Nw x 1)

% Collocation points (midpoints of each body panel)
co = 0.5*(Bp_w1 + Bp_w2);                        % Compute collocation point coordinates as (x,z) pairs
wake.co = co;                                % Collocation point coordinates (Nw x 2)
wake.cx = co(:,1);                           % Collocation point x coordinates (Nw x 1)
wake.cz = co(:,2);                           % Collocation point z coordinates (Nw x 1)
end