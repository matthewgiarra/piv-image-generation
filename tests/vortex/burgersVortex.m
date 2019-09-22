function [X, Y, Z, T] = burgersVortex(X0, Y0, Z0, TIMESPAN, PARAMS)
% [X, Y, Z, T] = burgersVortex(X0, Y0, TIMESPAN, PARAMS)
%
% PARAMS:
%   PARAMS.g = 10
%   PARAMS.a = 1;
%   PARAMS.m = 1e-6

% Constants for the burgers vortex
if isempty (PARAMS)
    PARAMS.g = 10;
    PARAMS.a = 1;
    PARAMS.m = 1e-6;
end

% Number of points
nPoints = numel(X0);

% Reshape the initial positions matrix
xo = reshape(X0, nPoints, 1);
yo = reshape(Y0, nPoints, 1);
zo = reshape(Z0, nPoints, 1);

% Solver options
solverOptions = odeset('RelTol', 1E-8, 'AbsTol', 1E-8);

% Run the solver
[t, positions] = ode45(@burgersVortexVelocityFunction, TIMESPAN,...
    [xo; yo; zo], solverOptions, PARAMS);

% If the specified time span was less than 2 elements long, take only the
% first and last elements of the solution. Otherwise return all of the elements.
% This is done because ODE45 calculates and returns many intermediate timesteps
% when the input time vector is only 2 elements long;
% If we didn't do this check here, then the returned solutions 
% would be at different time steps than were specified by the user in the JobFile.
if length(TIMESPAN) <= 2
    T = [t(1); t(end)];
    
    % Get [x,y,z] positions from positions vector
    x1 = positions(1,   1 : nPoints);
    x2 = positions(end, 1 : nPoints);
    
    y1 = positions(1,   nPoints + 1 : 2 * nPoints);
    y2 = positions(end, nPoints + 1 : 2 * nPoints);
    
    z1 = positions(1,   2 * nPoints + 1 : end);
    z2 = positions(end, 2 * nPoints + 1 : end);
    
    X = cat(1, x1, x2);
    Y = cat(1, y1, y2);
    Z = cat(1, z1, z2);
    
else
    T = t;
    X = positions(:, 1:nPoints);
    Y = positions(:, nPoints + 1 : 2*nPoints);
    Z = positions(:, 2*nPoints + 1 : end);
end

end



function U = burgersVortexVelocityFunction(T, X, PARAMS)

% X comes in as [x; y; z] where x, y, and z
% are column vectors of length nPoints
g = PARAMS.g;
a = PARAMS.a;
m = PARAMS.m;

nPoints = numel(X) / 3;
x = X(1 : nPoints);
y = X(nPoints + 1 : 2*nPoints);
z = X(2*nPoints + 1 : end);
[th, r, z] = cart2pol(x,y,z);

ur = -1/2 * a * r;
uth = g ./ (2 * pi * r) .* (1 - exp(-(a*r.^2) / (4 * m)));
uz = a * z;

ux = ur .* cos(th) - r .* uth .* sin(th);
uy = ur .* sin(th) + r .* uth .* cos(th);

U = [ux; uy; uz];        
        
end














