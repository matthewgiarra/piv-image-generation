function  [IMAGE1, IMAGE2] = generateImagePair_mc(...
    IMAGE_HEIGHT, IMAGE_WIDTH, PARTICLE_DIAMETER_MEAN, ...
    PARTICLE_DIAMETER_STD, PARTICLECONCENTRATION, ...
    DIFFUSION_STD_DEV, TRANSFORMATION_MATRIX, RUN_COMPILED)
% [IMAGE1, IMAGE2] = generateImagePair_mc(...
%     IMAGEHEIGHT, IMAGEWIDTH, PARTICLE_DIAMETER_MEAN, ...
%     PARTICLE_DIAMETER_STD, PARTICLECONCENTRATION, ...
%     DIFFUSION_STD_DEV, TRANSFORMATION_MATRIX, RUN_COMPILED)
% Generates a series of synthetic particle images that have been
%   transformed by a specified series of transformation matrices. 
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
%   DIFFUSION_STD_DEV = Standard deviation of the random component of each
%       particle's displacement, specified in pixels. The displacements
%       of the particles are calculated as the displacements due to the
%       specified affine transform plus a random displacement drawn
%       from a normal distribution with a mean of zero and a standard
%       deviation of DIFFUSION_STD_DEV.
%
%   TRANSFORMATIONMATRIX = 3 x 3 homogeneous matrix that specifies the
%   affine transformation relating the [x, y] positions of the particles
%   in the first image to the [x', y'] particle positions in the
%   second image. 
%
%   RUN_COMPILED = Boolean flag specifying whether or not to run the
%       compiled version of the particle image generation code.
%       The compiled executes about 100 times faster than the
%       non-compiled code.
%
% OUTPUTS
%   IMAGE_01 = 2D array containing the intensity values of the first
%   volumetric image.
%
%   IMAGE_02 = 2D array containing the intensity values of the second
%   volumetric image.
% 
% SEE ALSO
%    makeAffineTransform

% % % % % % % % % 
% BEGIN FUNCTION %
% % % % % % % % % 

% Double the height and width of the images so that rotations don't cause
% cropping.
augmentedWidth  = 2 * IMAGE_WIDTH; % Double the width so that rotations don't cause cropping
augmentedHeight = 2 * IMAGE_HEIGHT; % Double the height so that rotations don't cause cropping

% Calculate the coordinate of the geometric center of the image.
xc = (augmentedWidth + 1 ) / 2; % X center of image
yc = (augmentedHeight + 1 ) /2; % Y center of image

% Gaussian parameter to control the width of the
% simulated light-sheet intensity distribution in the depth direction
intensityStd = 0.25; % / (PARTICLEDIAMETER^2 / 8); 

% Number of particles to generate
nParticles = round(PARTICLECONCENTRATION * (augmentedHeight + 2 * PARTICLE_DIAMETER_MEAN) * (augmentedWidth + 2 * PARTICLE_DIAMETER_MEAN));

% Make random displacement matrix.
random_displacement_matrix = DIFFUSION_STD_DEV * randn([nParticles, 2]); 

% Uniformly distributed random numbers  from -0.5 to 0.5 corresponding to
% the depth-position of the particles in relation to the center of the 
% light sheet
particleCenters = rand( nParticles, 1 ) - 0.5; 

% Randomly generate the coordinates of particles
% Generate horizontal locations of particles in first image (column vector)
X1 = 1 + (augmentedWidth - 1) * rand(nParticles, 1);
Y1 = 1 + ( augmentedHeight - 1 ) * rand( nParticles, 1);

% Create a normal distribution of particle diameters
particle_diameters = abs(PARTICLE_DIAMETER_STD * randn(nParticles, 1) + PARTICLE_DIAMETER_MEAN);

% Gaussian Function that expresses the intensity distribution of the particle image on the "sensor"
particleMaxIntensities = exp( -particleCenters .^ 2 / (2 * intensityStd ^ 2 ) ); 

% Generate the first image. Choose between compiled and scripted code.
if RUN_COMPILED
    image1 = generateParticleImage_mex(augmentedHeight, augmentedWidth, X1, Y1, particle_diameters, particleMaxIntensities);
else
    image1 = generateParticleImage(    augmentedHeight, augmentedWidth, X1, Y1, particle_diameters, particleMaxIntensities);
end

% Crop the image and flip it vertically to place it in a right-handed coordinate system.
Image1Cropped = flipud(image1(augmentedHeight / 4 + 1 : 3 * augmentedHeight / 4, augmentedWidth / 4 + 1 : 3 * augmentedWidth / 4));

% Make 16 bit. Save to output variable.
IMAGE1 = uint16( (2^16 - 1) .* Image1Cropped .* 2.8^2 / PARTICLE_DIAMETER_MEAN ^2);

% Transform the particle coordinates
[Y2, X2] = transformImageCoordinates(TRANSFORMATION_MATRIX, X1, Y1, [yc xc]);

% Generate the second image. Choose between compiled and scripted code.
if RUN_COMPILED
    image2 = generateParticleImage_mex(augmentedHeight, augmentedWidth, X2 + random_displacement_matrix(:, 2), Y2 + random_displacement_matrix(:, 1), particle_diameters, particleMaxIntensities);
else
    image2 = generateParticleImage(    augmentedHeight, augmentedWidth, X2 + random_displacement_matrix(:, 2), Y2 + random_displacement_matrix(:, 1), particle_diameters, particleMaxIntensities);
end

% Crop the second image
Image2Cropped =  flipud(image2(augmentedHeight / 4 + 1: 3 * augmentedHeight / 4, augmentedWidth / 4 + 1: 3 * augmentedWidth / 4));

% Recast as 16-bit
IMAGE2 = uint16( (2^16 - 1) .* Image2Cropped .*  2.8^2 ./ PARTICLE_DIAMETER_MEAN^2);

end % End of function

