
function [X, Y, U, V, VORTICITY] = lambOseenVortexRingVorticityFunction(XIN, YIN, T, VORTEXPARAMETERS);

	if isempty(VORTEXPARAMETERS)
		
		
	end

circulation = VORTEXPARAMETERS.Circulation;
vortexRadius = VORTEXPARAMETERS.VortexRadius;
XC = VORTEXPARAMETERS.XC;
YC = VORTEXPARAMETERS.YC;
vortexAngle = VORTEXPARAMETERS.Angle;
viscosity = VORTEXPARAMETERS.Viscosity;


% % Convert input coordinates to column vectors
% x = X(1, :);
% y = X(2, :);

% Number of points
nPoints = numel(XIN);

x = XIN(:);
y = YIN(:);

% Positions of vortex cores (first core)
xc1 = XC - vortexRadius * sind(vortexAngle);
yc1 = YC + vortexRadius * cosd(vortexAngle);

% Positions of vortex cores (second core)
xc2 = XC + vortexRadius * sind(vortexAngle);
yc2 = YC - vortexRadius * cosd(vortexAngle);

% Transform the particle coordinates to polar coordinates
[th1, r1] = cart2pol(x - xc1 + 1, y - yc1 + 1);
[th2, r2] = cart2pol(x - xc2 + 1, y - yc2 + 1);

% Initialize vorticity
vorticity1 = zeros(nPoints, 1);
vorticity2 = zeros(nPoints, 1);

% % Initialize the angular positions of the particles in the second image
% uTheta1 = zeros(nPoints, 1);
% uTheta2 = zeros(nPoints, 1);

% Tangential velocities
uTheta1 =       circulation ./ (2 * pi * r1) .* (1 - exp(-r1.^2 ./ (4 * viscosity * T)));
uTheta2 = - 1 * circulation ./ (2 * pi * r2) .* (1 - exp(-r2.^2 ./ (4 * viscosity * T)));

% % Calculate tangential velocities for the first core
% vorticity1(r1 <= coreRadius) = circulation;
% vorticity1(r1 > coreRadius) =  0;
% 
% % Calculate tangential velocities for the second core
% vorticity2(r2 <= coreRadius) = -1 * circulation;
% vorticity2(r2 > coreRadius) =  0;

% Convert polar velocites to cartesian for first core
u1 = -1 * uTheta1 .* sin(th1);
v1 =      uTheta1 .* cos(th1);

% Convert polar velocites to cartesian for second core
u2 = -1 * uTheta2 .* sin(th2);
v2 =      uTheta2 .* cos(th2);

% Add the vorticies toget the total field.
VORTICITY = vorticity1 + vorticity2;

% Add the velocities together to get the total field. Reshape into row
% vectors.
U = reshape(u1 + u2, nPoints, 1);
V = reshape(v1 + v2, nPoints, 1);

% Name output variables.
X = x;
Y = y;

   
end