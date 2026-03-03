% Function to build panel mesh

% === INPUTS ===
% Bp - panel boundary points as (x,y) pairs ordered clockwise
% numPan - number of body panels

% === OUTPUTS ===
% == mesh - struct containing panel geometry data ==
% mesh.th - angle between panel and + x axis
% mesh.cos_th - cos(th) for each panel (useful for panel frame rotations) (numPan x 1)
% mesh.sin_th - sin(th) for each panel (useful for panel frame rotations) (numPan x 1)
% mesh.c - cos(th_i - th_j) or the cosine of the angle difference between each ith and jth panel (numPan x numPan)
% mesh.s - sin(th_i - th_j) or the cosine of the angle difference between each ith and jth panel (numPan x numPan)
% mesh.Bp1 - first NN-1 boundary points used for panel length and collocation point calculations (NN-1 x 2) 
% mesh.Bp2 - last NN-1 boundary points used for panel length and collocation point calculations (NN-1 x 2) 
% mesh.Bp1_x - x-coordinates of Bp1 (NN-1 x 1) 
% mesh.Bp1_z - z-coordinates of Bp1 (NN-1 x 1) 
% mesh.Bp2_x - x-coordinates of Bp2 (NN-1 x 1) 
% mesh.Bp2_z - z-coordinates of Bp2 (NN-1 x 1) 
% mesh.co - collocation point coordinates as (x,y) pairs (numPan x 2)
% mesh.cx - collocation point x-coordinates (numPan x 1)
% mesh.cz - collocation point z-coordinates (numPan x 1)

function mesh = build_panels(Bp, numPan)
% Build panel geometry, collocation points, and some reusable trig terms.

    % Panel vectors and panel angles (global frame)
    dxy     = diff(Bp);                          % (numPan x 2)
    th      = atan2(dxy(:,2), dxy(:,1));         % Panel angle wrt +x axis (numPan x 1)
    mesh.th = th;                                % Store panel angle
    mesh.cos_th = cos(th);                       % (numPan x 1) useful for rotating coordinate frame
    mesh.sin_th = sin(th);                       % (numPan x 1) useful for rotating coordinate frame

    % Precompute pairwise angle differences for rotating influence coefficients
    thj     = th(:).';                           % (1 x numPan)
    dth     = th - thj;                          % (numPan x numPan) via implicit expansion, angle differences between all panels
    mesh.c  = cos(dth);                          % cos(th_i - th_j)
    mesh.s  = sin(dth);                          % sin(th_i - th_j)

    % Panel endpoint vectors
    Bp1     = Bp(1:numPan, :);                  
    Bp2     = Bp(2:numPan+1, :);
    mesh.Bp = Bp;                                % Boundary point coordinates as (x,z) pairs (NN-1 x 2)
    mesh.Bp1 = Bp1;                              % Store Bp1
    mesh.Bp2 = Bp2;                              % Store Bp2
    mesh.Bp1_x = Bp1(:,1);  mesh.Bp1_z = Bp1(:,2); % Bp1 x and z coordinates (NN-1 x 1)
    mesh.Bp2_x = Bp2(:,1);  mesh.Bp2_z = Bp2(:,2); % Bp2 x and z coordinates (NN-1 x 1)

    % Collocation points (midpoints of each body panel)
    co = 0.5*(Bp1 + Bp2);                        % Compute collocation point coordinates as (x,z) pairs
    mesh.co = co;                                % Collocation point coordinates (numPan x 2)
    mesh.cx = co(:,1);                           % Collocation point x coordinates (numPan x 1)
    mesh.cz = co(:,2);                           % Collocation point z coordinates (numPan x 1)
   

end



