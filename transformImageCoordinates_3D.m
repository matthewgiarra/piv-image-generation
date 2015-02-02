%#codegen
function [YOUT, XOUT, ZOUT] = transformImageCoordinates_3D(TRANSFORM, XGRID, YGRID, ZGRID, CENTER)

% Count rows and columns in the original image 
[num_rows, num_cols, num_depths] = size(XGRID);

% Center of image
yc = CENTER(1);
xc = CENTER(2);
zc = CENTER(3);

% Shift grid origin to the center of grid
yPoints = YGRID - yc;
xPoints = XGRID - xc;
zPoints = ZGRID - zc;

% Reshape coordinate arrays into vectors
yPointsVect = reshape(yPoints, 1, numel(yPoints));
xPointsVect = reshape(xPoints, 1, numel(xPoints));
zPointsVect = reshape(zPoints, 1, numel(zPoints));

% Apply the transformation matrix to the shifted points
transformedPoints = TRANSFORM * [xPointsVect; yPointsVect; zPointsVect; ones( 1, length(xPointsVect ) ) ];

% Extract the X, Y, and Z coordinates of the transformed points.
Yt = transformedPoints(2, :);
Xt = transformedPoints(1, :);
Zt = transformedPoints(3, :);

% Turn transformed X, Y, and Z coordinates back into matrices and shift the origin back to (1, 1, 1)
YOUT = reshape(Yt, num_rows, num_cols, num_depths) + yc;
XOUT = reshape(Xt, num_rows, num_cols, num_depths) + xc;
ZOUT = reshape(Zt, num_rows, num_cols, num_depths) + yc;

end



