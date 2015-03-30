function [IMAGE_01, IMAGE_02] = generateImagePair_3D_mc(IMAGE_HEIGHT, ...
    IMAGE_WIDTH, IMAGE_DEPTH, PARTICLE_DIAMETER_MEAN, ...
    PARTICLE_DIAMETER_STD, PARTICLE_CONCENTRATION, BEAM_PLANE_STD_DEV, ...
    TRANSFORMATIONMATRIX, RUN_COMPILED)
%[IMAGE_01, IMAGE_02] = generateImagePair_3D_mc(IMAGE_HEIGHT, ...
%     IMAGE_WIDTH, IMAGE_DEPTH, PARTICLE_DIAMETER_MEAN, ...
%     PARTICLE_DIAMETER_STD, PARTICLE_CONCENTRATION, ...
%     TRANSFORMATIONMATRIX, RUN_COMPILED)
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
%   IMAGE_HEIGHT = Height of generated images in pixels
%
%   IMAGE_WIDTH = Width of generate images in pixels
%
%   IMAGE_DEPTH = Depth dimension of the images in pixels
%
%   PARTICLE_DIAMETER_MEAN = Mean of the particles' diameters in pixels.
%       Each particle's diameter will be drawn from a normal
%       distribution with a mean of PARTICLE_DIAMETER_MEAN and a standard
%       deviation of PARTICLE_DIAMETER_STD.
%
%   PARTICLE_DIAMETER_STD = Standard deviation of the particles' diameters
%       in pixels. Each particle's diameter will be drawn from a normal
%       distribution with a mean of PARTICLE_DIAMETER_MEAN and a standard
%       deviation of PARTICLE_DIAMETER_STD.
%
%   PARTICLE_CONCENTRATION = Image density of the particles in particles
%       per pixel.
%
%   TRANSFORMATIONMATRIX = 4 x 4 homogeneous matrix that specifies the
%   affine transformation relating the [x, y, z] positions of the particles
%   in the first image to the [x', y', z'] particle positions in the
%   second image. 
%
%   RUN_COMPILED = Boolean flag specifying whether or not to run the
%       compiled version of the particle image generation code.
%       The compiled executes about 100 times faster than the
%       non-compiled code.
%
% OUTPUTS
%   IMAGE_01 = 3D array containing the intensity values of the first
%   volumetric image.
%
%   IMAGE_02 = 3D array containing the intensity values of the second
%   volumetric image.
% 
% SEE ALSO
%    makeAffineTransform

% % % % % % % % % 
% BEGIN FUNCTION %
% % % % % % % % % 

% Double the height, width, and depth of the images so that rotations don't cause
% cropping.
augmentedWidth  = 2 * IMAGE_WIDTH;
augmentedHeight = 2 * IMAGE_HEIGHT;
augmentedDepth =  2 * IMAGE_DEPTH; 

% Calculate the coordinate of the geometric center of the image.
xc = (augmentedWidth  + 1 ) / 2; % X center of image
yc = (augmentedHeight + 1 ) / 2; % Y center of image
zc = (augmentedDepth  + 1 ) / 2; % Y center of image

% Number of particles to generate
nParticles = round(PARTICLE_CONCENTRATION * ...
    (augmentedHeight +  2 * PARTICLE_DIAMETER_MEAN) * ...
    (augmentedWidth + 2 * PARTICLE_DIAMETER_MEAN) * ...
    (augmentedDepth + 2 * PARTICLE_DIAMETER_MEAN));

% Randomly generate the coordinates of particles
% Generate horizontal locations of particles in first image (column vector)
X1 = 1 + (augmentedWidth  - 1) * rand(nParticles, 1);
Y1 = 1 + (augmentedHeight - 1) * rand(nParticles, 1);
Z1 = 1 + (augmentedDepth  - 1) * rand(nParticles, 1);

% Create a normal distribution of particle diameters
particle_diameters = abs(PARTICLE_DIAMETER_STD * randn(nParticles, 1) ...
    + PARTICLE_DIAMETER_MEAN);

% Gaussian Function that expresses the intensity distribution of the 
% particle image on the "sensor" This is set to 1 for now since the 
% depth-positions of the particles are known for 3D rendering.
% particleMaxIntensities = ones(nParticles, 1);

% Set the particle max intensities to be proportional to the 
% Gaussian intensity profile of the beam
particleMaxIntensities_01 = exp(-(Z1 - zc).^2 ./ ...
    (2 * BEAM_PLANE_STD_DEV ^ 2));

% Generate the first image.
% Choose whether to run compiled code.
if RUN_COMPILED
    image1 = generateParticleImage_3D_mex(augmentedHeight, ...
        augmentedWidth, augmentedDepth, X1, Y1, Z1, ...
        particle_diameters, particleMaxIntensities_01);
else
    image1 = generateParticleImage_3D(augmentedHeight, ...
        augmentedWidth, augmentedDepth, X1, Y1, Z1, ...
        particle_diameters, particleMaxIntensities_01);
end

% Crop the image and flip it vertically to place it in a right-handed coordinate system.
Image1Cropped = flipud(image1(augmentedHeight   / 4 + ...
    1 : 3 * augmentedHeight / 4, augmentedWidth / 4 + ...
    1 : 3 * augmentedWidth  / 4, augmentedDepth / 4 + ...
    1 : 3 * augmentedDepth  / 4));

% Make 16 bit. Save to output variable.
IMAGE_01 = uint16( (2^16 - 1) .* ...
    Image1Cropped .* 2.8 ^ 2 / PARTICLE_DIAMETER_MEAN ^ 2);

% Transform the particle coordinates
[Y2, X2, Z2] = transformImageCoordinates_3D(TRANSFORMATIONMATRIX, ...
    X1, Y1, Z1, [yc, xc, zc]);

% Set the particle max intensities to be proportional to the 
% Gaussian intensity profile of the beam
particleMaxIntensities_02 = exp(-(Z2 - zc).^2 / ...
    (2 * BEAM_PLANE_STD_DEV ^ 2));


% Generate the second image
% Choose whether to run compiled code.
if RUN_COMPILED
    image2 = generateParticleImage_3D_mex(augmentedHeight, ...
        augmentedWidth, augmentedDepth, X2, Y2, Z2, particle_diameters, ...
        particleMaxIntensities_02);
else
    image2 = generateParticleImage_3D(augmentedHeight, ...
        augmentedWidth, augmentedDepth, X2, Y2, Z2, particle_diameters, ...
        particleMaxIntensities_02);
end

% Crop the second image
Image2Cropped =  flipud(image2(augmentedHeight  / 4 + ...
    1 : 3 * augmentedHeight / 4, augmentedWidth / 4 + ...
    1 : 3 * augmentedWidth  / 4, augmentedDepth / 4 + ...
    1 : 3 * augmentedDepth  / 4));

% Recast as 16-bit
IMAGE_02 = uint16( (2^16 - 1) .* Image2Cropped .* ...
    2.8 ^ 2 ./ PARTICLE_DIAMETER_MEAN ^ 2);

end % End of function

