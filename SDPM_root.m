% Wrapper script to execute Dirichlet_ConstantSourceDoublet.m over a
% specified number of iterations (time-steps)

UserInput.steady = 'steady'
UserInput.timeStep = 0.1 % Time increments [s] at which the Panel code is executed 
UserInput.Freestream = 1; % Free-stream velocity magnitude
UserInput.alpha = 5; % Angle of attack [degrees] (angle of attack at t = 0 for unsteady case)
