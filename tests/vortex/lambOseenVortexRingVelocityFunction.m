function [U, VORTICITY] = lambOseenVortexRingVelocityFunction(T, X, VORTEXPARAMETERS);
% [U, VORTICITY] = lambOseenVortexRingVelocityFunction(T, X, VORTEXPARAMETERS)
% This function calculates the analytical velocity and vorticity of a Lamb-Oseen vortex ring
% at specified positions. Note that this solution does not include diffusion of momentum due to viscosity; 
% This means that the vortex ring does not decay over time. Also note that
% the coordinate system is right-handed; i.e., the point [1,1] corresponds
% to the lower left corner of the field. 
%
% INPUTS
% T = Time at which to evaluate the solution (in frames). When calling this
%   function independently (i.e., when it is not called by ODE45), T should
%   be a single element, not a vector. When called by ODE45, T should be a
%   column vector of times.
%
% X = Column vector of horizontal and vertical coordinates at which to 
%   evaluate the solution (in pixels). 
%   This formatting is required by ODE45. The length of the vector should
%   be twice the number of coordinates (nPoints); The horizontal positions
%   are specified by the elements X(1 : nPoints) and the corresponding
%   vertical positions are specified by the elements X(nPoints + 1 : end).
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
%   VORTEXPARAMETERS.XC = Horizontal position of the center of the vortex ring (pixels)
% 
%   VORTEXPARAMETERS.YC = Vertical position of the center of the vortex ring (pixels)
%
% OUTPUTS
%   U = Column vector specifying the horizontal and vertical velocities (pixels/frame)
%   evaluated at time T.  
%   This formatting is required by ODE45. The length of the vector should
%   be twice the number of coordinates (nPoints); The horizontal velocities
%   are specified by the elements U(1 : nPoints) and the corresponding
%   vertical velocities are specified by the elements U(nPoints + 1 : end).
%
%   VORTICITY = Column vector specifying the vorticity at each coordinate
%   (in 1/frames). 
%
% EXAMPLE
% T = 0;
% [x, y] = meshgrid(1:1024, 1:1024);
% nPoints = numel(x);
% X = cat(1, x(:), y(:));
% VORTEXPARAMETERS.CoreRadius = 100;
% VORTEXPARAMETERS.VortexRadius = 200;
% VORTEXPARAMETERS.Angle = 0;
% VORTEXPARAMETERS.PeakVelocity = 5;
% VORTEXPARAMETERS.PropagationVelocity = 0;
% VORTEXPARAMETERS.XC = 512;
% VORTEXPARAMETERS.YC = 512;
%
% % Calculate the velocity and vorticity fields
% [U, VORTICITY] = lambOseenVortexRingVelocityFunction(T, X, VORTEXPARAMETERS);
%
% % Extract the velocities and reshape into matrices.
% u = reshape(U(1:nPoints), 1024, 1024);
% v = reshape(U(nPoints + 1 : end), 1024, 1024);
% w = reshape(VORTICITY, 1024, 1024);
%
% % Plot the velocity and vorticity fields.
% contourf(x, y, w);
% hold on;
% quiver(x(1:10:end, 1:10:end), y(1:10:end, 1:10:end), u(1:10:end, 1:10:end), v(1:10:end, 1:10:end), 1.5, 'black');
% colorbar;
% t = colorbar('peer', gca);
% set(get(t, 'ylabel'), 'String', 'Vorticity (1/frames)', 'FontSize', 18);
% axis image;
% set(gcf, 'color', 'white');
% set(gca, 'FontSize', 18);
% title('Lamb-Oseen Vortex Ring (velocity and vorticity)');
% hold off

% Begin Function

% Count the number of points
nPoints = numel(X)/2;

% Extract the horizontal and vertical coordinates from the column-vector
% input. The input is a column vector because it has to work with ODE45,
% which demands a column vector input.
x = X(1 : nPoints);
y = X(nPoints + 1 : end);

% Check to see if vortex parameters were supplied
if nargin == 2
	VORTEXPARAMETERS = [];
end

% Default options
if isempty(VORTEXPARAMETERS)
	
	% Min and max x
	x_min = min(x(:));
	x_max = max(x(:));
	
	% Min and max y
	y_min = min(y(:));
	y_max = max(y(:));
	
	% Guess the height and width
	domain_height = y_max  +  (y_min < 0) * y_min;
	domain_width  = x_max  +  (x_min < 0) * x_min;
	
	% Minimum dimension
	min_dim = min([domain_height, domain_width]);
	
	VORTEXPARAMETERS.CoreRadius = min_dim / 5;
	VORTEXPARAMETERS.VortexRadius = domain_height / 5;
	VORTEXPARAMETERS.Angle = 0;
	VORTEXPARAMETERS.PeakVelocity = 5;
	VORTEXPARAMETERS.PropagationVelocity = 0;
	VORTEXPARAMETERS.XC = domain_width / 2;
	VORTEXPARAMETERS.YC = domain_height / 2;
end

% Extract vortex parameters from structure
vortexRadius = VORTEXPARAMETERS.VortexRadius;
coreRadius = VORTEXPARAMETERS.CoreRadius;
XC = VORTEXPARAMETERS.XC;
YC = VORTEXPARAMETERS.YC;
vortexAngle = VORTEXPARAMETERS.Angle;
propagationVelocity = VORTEXPARAMETERS.PropagationVelocity;
maxTangentialVelocity = VORTEXPARAMETERS.PeakVelocity;



% Center of the vortex ring propagated by the vortex velocitity
xc = XC + propagationVelocity * T * cosd(vortexAngle);
yc = YC + propagationVelocity * T * sind(vortexAngle);

% Positions of vortex cores (first core)
xc1 = xc - vortexRadius * sind(vortexAngle);
yc1 = yc + vortexRadius * cosd(vortexAngle);

% Positions of vortex cores (second core)
xc2 = xc + vortexRadius * sind(vortexAngle);
yc2 = yc - vortexRadius * cosd(vortexAngle);

% % Transform the particle coordinates to polar coordinates
[th1, r1] = cart2pol(x - xc1 + 1, y - yc1 + 1);
[th2, r2] = cart2pol(x - xc2 + 1, y - yc2 + 1);

% Transform the particle coordinates to polar coordinates
% [th1, r1] = cart2pol(x - xc1, y - yc1);
% [th2, r2] = cart2pol(x - xc2, y - yc2);

% Parameter ("alpha") defined by Davenport et al, 1996, "The
% structure and development of a wing-tip vortex"
alpha = 1.25643;

% Calculate the tangential velocities due to each vortex core. 
uTheta1 =      maxTangentialVelocity * (1 + 0.5 / alpha ) * (coreRadius ./ r1) .* (1 - exp(-1 * alpha  .* r1.^2 / coreRadius^2));
uTheta2 = -1 * maxTangentialVelocity * (1 + 0.5 / alpha ) * (coreRadius ./ r2) .* (1 - exp(-1 * alpha  .* r2.^2 / coreRadius^2));

% Set NaN velocities (which occur at zero radius) equal to zero
uTheta1(isnan(uTheta1)) = 0;
uTheta2(isnan(uTheta2)) = 0;

% Calculate vorticity for the first core
vort1 =      maxTangentialVelocity * coreRadius * (1 + 0.5 / alpha) * 2 * alpha / coreRadius^2 .* exp(-1 * alpha .* r1.^2 / coreRadius^2);
vort2 = -1 * maxTangentialVelocity * coreRadius * (1 + 0.5 / alpha) * 2 * alpha / coreRadius^2 .* exp(-1 * alpha .* r2.^2 / coreRadius^2);

% Add the vorticities
VORTICITY = vort1 + vort2;

% Convert polar velocites to cartesian for first core
u1 = -1 * uTheta1 .* sin(th1);
v1 =      uTheta1 .* cos(th1);

% Convert polar velocites to cartesian for second core
u2 = -1 * uTheta2 .* sin(th2);
v2 =      uTheta2 .* cos(th2);

% Add the velocities together to get the total field. Reshape into column vectors.
u = reshape(u1 + u2, nPoints, 1);
v = reshape(v1 + v2, nPoints, 1);

% Concatonate into a single vector.
U = [u; v];
   
end