%#codegen
function [X_OUT, Y_OUT, Z_OUT] = transformImageCoordinates_3D(...
    X_IN, Y_IN, Z_IN, TRANSFORM)

% Measure size of the input coordinates
[num_rows, num_cols, depth] = size(X_IN);

% Apply the transformation matrix to the shifted points
transformedPoints = TRANSFORM * [(X_IN(:))'; (Y_IN(:))'; (Z_IN(:))'; ...
    ones( 1, numel(X_IN)) ];

% Transformed X coordinate
X_OUT = reshape(transformedPoints(1, :), ...
    [num_rows, num_cols, depth]);

% Transformed Y coordinate
Y_OUT = reshape(transformedPoints(2, :), ...
    num_rows, num_cols, depth);

% Transformed Z coordinate
Z_OUT = reshape(transformedPoints(3, :), ...
    num_rows, num_cols, depth);

end



