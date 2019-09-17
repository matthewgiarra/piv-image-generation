
function U = burgersVortexVelocityFunction(T, X, PARAMS);

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

% f = @(t,x) [(-1/2)*A*r; 
%             (G/(2*pi*r))*(1-exp(-(A*r^2)/(4*m)));
%             A*z];
        
        
        
end