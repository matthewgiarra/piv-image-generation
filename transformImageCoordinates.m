%#codegen
function [YOUT, XOUT] = transformImageCoordinates(TRANSFORM, XGRID, YGRID, CENTER)

% Count rows and columns in the original image 
[nr, nc] = size(XGRID);

% Center of image
yc = CENTER(1);
xc = CENTER(2);

% Shift grid origin to the center of grid
% Is this affecting things???
xPoints = XGRID - xc;
yPoints = YGRID - yc;

% Reshape coordinate matrices into vectors
xPointsVect = reshape(xPoints, 1, numel(xPoints));
yPointsVect = reshape(yPoints, 1, numel(yPoints));

% Apply the transformation matrix to the shifted points
transformedPoints = TRANSFORM * [xPointsVect; yPointsVect; ones( 1, length(xPointsVect ) ) ];

% Extract the X- and Y coordinates of the transformed points.
Xt = transformedPoints(1, :);
Yt = transformedPoints(2, :);

% Turn transformed X- and Y- coordinates back into matrices and shift the origin back to (1, 1)
% x = reshape(Xt, nr, nc) + xc;
% y = reshape(Yt, nr, nc) + yc;
XOUT = reshape(Xt, nr, nc) + xc;
YOUT = reshape(Yt, nr, nc) + yc;


% % Crop the image if specified.
% if strcmp(METHOD, 'crop')
%     XOUT = (x(x >= min(XGRID(:)) & x <= max(XGRID(:)) & y >= min(YGRID(:)) & y <= max(YGRID(:))));
%     YOUT = (y(x >= min(XGRID(:)) & x <= max(XGRID(:)) & y >= min(YGRID(:)) & y <= max(YGRID(:))));
%     
% else % Don't crop the image if cropping wasn't specified.
%     XOUT = x;
%     YOUT = y;
%     
% end


end



