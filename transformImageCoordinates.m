%#codegen
function [YOUT, XOUT] = transformImageCoordinates(TRANSFORM, XGRID, YGRID, CENTER)

% Count rows and columns in the original image 
[num_rows, num_cols] = size(XGRID);

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
XOUT = reshape(Xt, num_rows, num_cols) + xc;
YOUT = reshape(Yt, num_rows, num_cols) + yc;

end



