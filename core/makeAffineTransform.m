function [TRANSFORM, SCALING, ROTATION, SHEARING, TRANSLATION] = ...
    makeAffineTransform(SCALE, ROTATIONANGLE, SHEARX, SHEARY, ...
    TRANSLATIONX, TRANSLATIONY, PLOT)
% MAKEAFFINETRANSFORM(SCALEX, SCALEY, ROTATIONANGLE, SHEARX, SHEARY, TRANSLATIONX, TRANSLATIONY, PLOT)
% Creates a matrix that performs an affine transformation on a vector or set of vectors; that is, it performs scaling,
% rotation, shearing, and translation. 
%
% INPUTS
%   SCALE = Uniform Scaling factor
%   ROTATIONANGLE = Rotation angle from the positive horizontal axis in a
%       right-handed coordinate system (radians). A positive rotation angle specifies a
%       counter-clockwise rotation.
%   SHEARX = Horizontal shear rate
%   SHEARY = Vertical shear rate
%   TRANSLATIONX = Horizontal translation
%   TRANSLATIONY = Vertical translation
%   PLOT = Flag to specify whether or not to produce a figure to
%       illustrate the effect of the transformation matrix on a square
%       quadrilateral (1 to generate plot; 0 to suppress plot).
% 
% OUTPUTS
%   TRANSFORM = 3 x 3 homogeneous matrix specifying the total affine transformation.
%   SCALING = 3 x 3 matrix specifying the pure scaling transformation.
%   ROTATION = 3 x 3 matrix specifying the pure rotation transformation. 
%   SHEARING = 3 x 3 matrix specifying the pure shearing transformation.
%   TRANSLATION = 3 x 3 matrix specifying the pure translation transformation.
% 

% %%%%%%%%%%
% BEGIN FUNCTION %
%%%%%%%%%%%

% Default to not generate plots
if nargin < 7 
    PLOT = 0;
end

% Isotropic scaling matrix
SCALING = [SCALE 0 0; 0 SCALE 0; 0 0 1]; 

% Rotation matrix
ROTATION = [cos(ROTATIONANGLE) -sin(ROTATIONANGLE) 0; sin(ROTATIONANGLE) cos(ROTATIONANGLE) 0; 0 0 1]; 

% Shearing matrix
SHEARING = [1 SHEARX 0; SHEARY 1 0; 0 0 1]; 

% Translation matrix
TRANSLATION = [1 0 TRANSLATIONX; 0 1 TRANSLATIONY; 0 0 1]; 

% Affine transformation matrix 
TRANSFORM = ( TRANSLATION * SCALING * ROTATION * SHEARING ); 

if PLOT % If the "plot" flag is enabled, generate a figure to illustrate the effect of the transformation matrix on a square.
%     The positive direction of the vertical (Y) axis is downward for consistency with image axes.
    X1 = [-1 1 1 -1]'; % X-coordinates of a test quadrilateral 
    Y1 = [-1 -1 1 1]'; % Y-coordinates of a test quadrilateral 
    Z2 = [X1 Y1 ones(length(X1), 1)] * TRANSFORM ; % Apply the transformation to the undeformed coordinates (i.e., deform them).
    X2 = Z2(:, 1); % Deformed X-coordinates
    Y2 = Z2(:, 2); % Deformed Y-coordinates
    patch(X1, Y1, 'green'); hold on; patch(X2, Y2, 'blue', 'FaceAlpha', 0.5); % Overlay the deformed and undeformed polygons
    axis image; set(gcf, 'Color', [1 1 1]); set(gca, 'FontSize', 16) % Format the plot
    set(gca, 'YDir', 'reverse');
    leg = legend('Undeformed', 'Deformed', 'location', 'NorthEastOutside'); set(leg, 'FontSize', 16); % Format the legend
    hold off; % Release the figure hold.
end

end

% %%%%%%%%%%
% END OF FUNCTION %
%%%%%%%%%%%
