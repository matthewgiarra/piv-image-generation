function [TRANSFORM, SCALING, ROTATION_ANGLES, TRANSLATION] = ...
    makeSimilarityTransform_3D(SCALE, ROTATION_ANGLES, ...
    TRANSLATION_X, TRANSLATION_Y, TRANSLATION_Z)
% makeSimilarityTransform_3D(SCALE ROTATION_ANGLES, TRANSLATION_X, TRANSLATION_Y, TRANSLATION_Z);
% Creates a matrix that performs an affine transformation on a vector or set of vectors; that is, it performs scaling,
% rotation, shearing, and translation. 
%
% NOTE: This function needs to be checked against the angle measurement
% system used in the SOFT code.
%
% INPUTS
%   SCALE = Uniform Scaling factor
%   ROTATION_ANGLES = Three-element vector containing the Euler angles 
%       that describe the 3D rotation. See Kostelec 2008. 
%   TRANSLATION_X = Horizontal translation
%   TRANSLATION_Y = Vertical translation
%   TRANSLATION_Z = Translation in the Z direction
% 
% OUTPUTS
%   TRANSFORM = 4 x 4 homogeneous matrix specifying the 3D similarity transformation.
%   This matrix is compatible with the Matlab function IMTRANSFORM (i.e.,
%       the translation elements are in the third row rather than in the third column). 
%   SCALING = 3 x 3 matrix specifying the pure scaling transformation.
%   ROTATION = 3 x 3 matrix specifying the pure rotation transformation. 
%   SHEARING = 3 x 3 matrix specifying the pure shearing transformation.
%   TRANSLATION = 3 x 3 matrix specifying the pure translation transformation.
% 

% %%%%%%%%%%
% BEGIN FUNCTION %
%%%%%%%%%%%

% Isotropic scaling matrix (homogeneous form)
SCALING = [SCALE, 0, 0, 0; ...
           0, SCALE, 0, 0; ...
           0, 0, SCALE, 0; ...
           0, 0, 0, 1]; 

% This is the rotation due to the first Euler angle, which 
% is rotation about the Z axis
R_01 = [cos(ROTATION_ANGLES(1)), -sin(ROTATION_ANGLES(1)), 0, 0; ...
               sin(ROTATION_ANGLES(1)),  cos(ROTATION_ANGLES(1)), 0, 0; ...
               0, 0, 1, 0; ...
               0, 0, 0, 1];

% This is the rotation due to the second Euler angle, which 
% is rotation about the y axis
R_02 = [cos(ROTATION_ANGLES(2)), 0, sin(ROTATION_ANGLES(2)), 0; ...
                0, 1, 0, 0; ...
              -sin(ROTATION_ANGLES(2)), 0, cos(ROTATION_ANGLES(2)), 0; ...
                0, 0, 0, 1];
            
% This is the rotation due to the third Euler angle, which 
% is rotation about the new Z axis (I think)            
R_03 = [cos(ROTATION_ANGLES(3)), -sin(ROTATION_ANGLES(3)), 0, 0; ...
               sin(ROTATION_ANGLES(3)),  cos(ROTATION_ANGLES(3)), 0, 0; ...
               0, 0, 1, 0; ...
               0, 0, 0, 1];
 
% Combine the decomposed rotation matrices into
% a single homogeneous rotation matrix.
ROTATION = R_01 * R_02 * R_03;
           
% Translation matrix
TRANSLATION = [1, 0, 0, TRANSLATION_X; ...
               0, 1, 0, TRANSLATION_Y; ...
               0, 0, 1, TRANSLATION_Z; ...
               0, 0, 0, 1]; 

% Affine transformation matrix 
TRANSFORM = (TRANSLATION * SCALING * ROTATION); 

end

% %%%%%%%%%%
% END OF FUNCTION %
%%%%%%%%%%%
