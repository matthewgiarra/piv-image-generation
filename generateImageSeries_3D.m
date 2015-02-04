%#codegen
function IMAGE_SERIES = generateImageSeries_3D(IMAGE_HEIGHT, IMAGE_WIDTH, IMAGE_DEPTH, NUM_IMAGES, PARTICLE_DIAMETER_MEAN, PARTICLE_DIAMETER_STD, PARTICLECONCENTRATION, TRANSFORMS, COMPILED)
% generateImageSeries_3D(IMAGE_HEIGHT, IMAGE_WIDTH, IMAGE_DEPTH, NUM_IMAGES, PARTICLE_DIAMETER_MEAN, PARTICLE_DIAMETER_STD, PARTICLECONCENTRATION, TRANSFORMS, COMPILED)
% Generates a series of synthetic particle images that have been transformed by a
% specified series of transformation matrices. 
%
% The equations used to generate the particles were taken from Equation 16 (on page 5) of the paper
% "Methods for Digital Particle Image Sizing (DPIS): Comparisons and
% Improvements" (2009) by Michael R. Brady, Samuel G. Raben, Pavlos P.
% Vlachos, published in the journal Flow Measurement and Instrumentation.
% (doi:10.1016/j.flowmeasinst.2009.08.001)
%
% INPUTS
%   IMAGE_HEIGHT = Height of generated images (pixels)
%
%   IMAGE_WIDTH = Width of generate images (pixels)
%
%   PARTICLE_CONCENTRATION = Average number of particles per pixel
%
%   TRANSFORMATIONMATRIX = 3 x 3 x NUM_IMAGES matrix that specifies the affine
%   transformation relating the particle positions in the first image to
%   the particle positions in the second image. 
%
%   SAVEFLAG = Flag to specify whether or not to save the generated images (1 to save, 0 to not save). 
%
%   SAVEDIR = Directory in which to save the generated images. 
%
%   SAVENAME = Name of the images to save, which will be followed by '_A'
%   and '_B' for the first and second images in the pair. 
%
%   CONTROLCASE = Binary flag (0 or 1) specifying whether to generate a
%   simple six-particle pattern to easily observe the effects of the
%   transformation on the images. This should only be enabled for debugging
%   purposes. 
%
% OUTPUTS
%   This code does not save any outputs to the workspace. Instead, it saves
%   images to locations specified by the inputs to the code.
% 
% SEE ALSO
%    makeAffineTransform

% % % % % % % % % 
% BEGIN FUNCTION %
% % % % % % % % % 

% Default to not running compiled codes.
if nargin < 9
    COMPILED = 0;
end;

% Allocate the image list
IMAGE_SERIES = zeros(IMAGE_HEIGHT, IMAGE_WIDTH, IMAGE_DEPTH, NUM_IMAGES, 'uint16');

% Double the height, width, and depth of the images so that rotations don't cause
% cropping.
augmentedWidth  = 2 * IMAGE_WIDTH;
augmentedHeight = 2 * IMAGE_HEIGHT;
augmentedDepth =  2 * IMAGE_DEPTH; 

% Calculate the coordinate of the geometric center of the image.
xc = (augmentedWidth + 1 )  / 2; % X center of image
yc = (augmentedHeight + 1 ) / 2; % Y center of image
zc = (augmentedDepth + 1 )  / 2; % Y center of image

% Number of particles to generate
nParticles = round(PARTICLECONCENTRATION * (augmentedHeight + 2 * PARTICLE_DIAMETER_MEAN) * (augmentedWidth + 2 * PARTICLE_DIAMETER_MEAN) * (augmentedDepth + 2 * PARTICLE_DIAMETER_MEAN));

% Randomly generate the coordinates of particles
% Generate horizontal locations of particles in first image (column vector)
X1 = 1 + (augmentedWidth  - 1) * rand(nParticles, 1);
Y1 = 1 + (augmentedHeight - 1) * rand(nParticles, 1);
Z1 = 1 + (augmentedDepth  - 1) * rand(nParticles, 1);

% Create a normal distribution of particle diameters
particle_diameters = abs(PARTICLE_DIAMETER_STD * randn(nParticles, 1) + PARTICLE_DIAMETER_MEAN);

% Gaussian Function that expresses the intensity distribution of the particle image on the "sensor"
% This is set to 1 for now since the depth-positions of the particles are
% known for 3D rendering.
particleMaxIntensities = ones(nParticles, 1);

% Create a placeholder to store the generated image
if COMPILED
    % Run compiled MEX code
    image_01 = generateParticleImage_3D_mex(augmentedHeight, augmentedWidth, augmentedDepth, X1, Y1, Z1, particle_diameters, particleMaxIntensities);
else
    % Run scripted Matlab code
    image_01 = generateParticleImage_3D(augmentedHeight, augmentedWidth, augmentedDepth, X1, Y1, Z1, particle_diameters, particleMaxIntensities);
end

% Crop the image and flip it vertically to place it in a right-handed coordinate system.
Image1Cropped = flipud(image_01(augmentedHeight / 4 + 1 : 3 * augmentedHeight / 4, augmentedWidth / 4 + 1 : 3 * augmentedWidth / 4, augmentedDepth / 4 + 1 : 3 * augmentedDepth / 4));

% Make 16 bit. Save to output variable.
IMAGE_SERIES(:, :, :, 1) = uint16( (2^16 - 1) .* Image1Cropped .* 2.8^2 / PARTICLE_DIAMETER_MEAN ^2);

% Generate the rest of the images, which are
% images containing particles whose coordinates
% are transformed copies of the coordinates of the
% particles in the first image.
for k = 2 : NUM_IMAGES
    
    % Transform the particle coordinates
    [Y2, X2, Z2] = transformImageCoordinates_3D(TRANSFORMS(:, :, k), X1, Y1, Z1, [yc, xc, zc]);

    % Generate the second image
    if COMPILED
        % Run compiled MEX code
        image_02 = generateParticleImage_3D_mex(augmentedHeight, augmentedWidth, augmentedDepth, X2, Y2, Z2, particle_diameters, particleMaxIntensities);
    else
        % Generate the second image
        image_02 = generateParticleImage_3D(augmentedHeight, augmentedWidth, augmentedDepth, X2, Y2, Z2, particle_diameters, particleMaxIntensities);
    end
   
    % Crop the second image
    Image2Cropped =  flipud(image_02(augmentedHeight / 4 + 1: 3 * augmentedHeight / 4, augmentedWidth / 4 + 1: 3 * augmentedWidth / 4, augmentedDepth / 4 + 1: 3 * augmentedDepth / 4));

    % Recast as 16-bit
    IMAGE_SERIES(:, :, :, k) = uint16( (2^16 - 1) .* Image2Cropped .*  2.8^2 ./ PARTICLE_DIAMETER_MEAN^2);

end

end % End of function

