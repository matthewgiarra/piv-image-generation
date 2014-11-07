function [X, Y, T] = lambOseenVortexRing(X0, Y0, TIMESPAN, VORTEXPARAMETERS);
% [X, Y, T] = lambOseenVortexRing(X0, Y0, TIMESPAN, VORTEXPARAMETERS)
% This function advects massless flow tracers with initial positions [X0, Y0]
% over a time-vector specified by TIMESPAN according to a Lamb-Oseen vortex
% ring velocity field whose parameters are specified by VORTEXPARAMETERS. 
% Note that this solution does not include diffusion of momentum due to viscosity; 
% This means that the vortex ring does not decay over time. The positions
% of flow tracers are integrated using ODE45. 
%
% INPUTS
% X0 = Matrix or vector of initial horizontal positions flow tracers
%   (in pixels; can be non-integer)
% Y0 = Matrix or vector of initial vertical positions flow tracers
%   (in pixels; can be non-integer)
% TIMESPAN = Column vector of times at which to evaluate the particle
%   positions (in frames)
%
% VORTEXPARAMETERS
%   This is a structure whose fields specify the parameters of the vortex
%       ring. Each field is described below.
%
%   VORTEXPARAMETERS.CoreRadius = Distance from the center of each vortex
%       core to the max tangential velocity produced by that core (pixels)
%
%   VORTEXPARAMETERS.VortexRadius = Distance between the centers 
%       of the two vortex cores (pixels);
%
%   VORTEXPARAMETERS.Angle = Angle (in degrees) between the direction of propagation and
%       the positive horizontal axis; positive angles indicate counterclockwise
%       rotation, and 0 angle indicates propagation to the right. 
%
%   VORTEXPARAMETERS.PeakVelocity = Peak velocity produced by each vortex core (pixels/frame);
%
%   VORTEXPARAMETERS.PropagationVelocity = Propagation velocity of the vortex ring (pixels/frame);
%
% OUTPUTS
%   X = Matrix of horizontal positions of flow tracers at each solution
%   time. Each row of this matrix corresponds to a single time instance of
%   the solution; the k'th row corresponds to the solution evaluated at the
%   time specified by the k'th element of the input variable TIMESPAN.
%
%   Y = Matrix of vertical positions of flow tracers at each solution
%   time. Each row of this matrix corresponds to a single time instance of
%   the solution; the k'th row corresponds to the solution evaluated at the
%   time specified by the k'th element of the input variable TIMESPAN.
%
%   T = Times at which the solution is evaluated; this vector is identical
%   to the input TIMESPAN.
%
% EXAMPLE
% % Generate random horizontal positions between 1 and 1024;
% X0 = 1 + 1023 * rand(1000, 1);
% 
% % Generate random vertical positions between 1 and 1024;
% Y0 = 1 + 1023 * rand(1000, 1);
% 
% % Time vector
% TIMESPAN = [0, 2];
% 
% % Set vortex parameters
% VORTEXPARAMETERS.CoreRadius = 100;
% VORTEXPARAMETERS.VortexRadius = 200;
% VORTEXPARAMETERS.Angle = 0;
% VORTEXPARAMETERS.PeakVelocity = 5;
% VORTEXPARAMETERS.PropagationVelocity = 0;
% VORTEXPARAMETERS.XC = 512;
% VORTEXPARAMETERS.YC = 512;
% 
% % Integrate particle positions
% [X, Y, T] = lambOseenVortexRing(X0, Y0, TIMESPAN, VORTEXPARAMETERS);
% 
% % Make a plot of initial and final particle positions
% figure(1)
% plot(X(1, :), Y(1, :), 'ok', 'MarkerFaceColor', 'black'); axis image
% hold on
% plot(X(end, :), Y(end, :), 'or', 'MarkerFaceColor', 'red');
% title('Particle positions', 'FontSize', 18);
% h = legend('Initial Positions', 'Final Positions');
% axis image
% set(h, 'FontSize', 18);
% set(gcf, 'Color', 'white');
% set(gca, 'FontSize', 18);
% hold off;
% 
% % Make a plot of the displacements
% figure(2);
% U = X(end, :) - X(1, :);
% V = Y(end, :) - Y(1, :);
% quiver(X(1, :), Y(1, :), U, V, 'black');
% title('Displacement vectors', 'FontSize', 18);
% axis image;
% set(gcf, 'color', 'white');
% set(gca, 'FontSize', 18);

% Number of points
nPoints = numel(X0);

% Reshape the initial positions matrix
xo = reshape(X0, nPoints, 1);
yo = reshape(Y0, nPoints, 1);

% Solver options
% solverOptions = [];
solverOptions = odeset('RelTol', 1E-8, 'AbsTol', 1E-8);

% Run the solver
[t, positions] = ode45(@lambOseenVortexRingVelocityFunction, TIMESPAN, [xo; yo], solverOptions, VORTEXPARAMETERS);

% If the specified time span was less than 2 elements long, take only the
% first and last elements of the solution. Otherwise return all of the elements.
% This is done because ODE45 calculates and returns many intermediate timesteps
% when the input time vector is only 2 elements long;
% If we didn't do this check here, then the returned solutions 
% would be at different time steps than were specified by the user in the JobFile.
if length(TIMESPAN) <= 2
    T = [t(1); t(end)];
    X = cat(1, positions(1, 1:nPoints), positions(end, 1:nPoints));
    Y = cat(1, positions(1, nPoints + 1 : end), positions(end, nPoints + 1 : end));
else
    T = t;
    X = positions(:, 1:nPoints);
    Y = positions(:, nPoints + 1 : end);
end

end















