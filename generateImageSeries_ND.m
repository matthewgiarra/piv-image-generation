function IMAGE_SERIES = generateImageSeries_ND(IMAGE_HEIGHT, ...
    IMAGE_WIDTH, IMAGE_DEPTH, PARTICLE_DIAMETER_MEAN, ...
    PARTICLE_DIAMETER_STD, PARTICLECONCENTRATION, ...
    BEAM_PLANE_STD_DEV, TRANSFORMS, DIMS, PARALLEL_PROCESSING, COMPILED)
% generateImageSeries_3D(IMAGE_HEIGHT, IMAGE_WIDTH, 
% IMAGE_DEPTH, NUM_IMAGES, PARTICLE_DIAMETER_MEAN, PARTICLE_DIAMETER_STD, 
% PARTICLECONCENTRATION, TRANSFORMS, DIMS, COMPILED)
% Generates a series of synthetic particle images that have
% been transformed by a specified series of transformation matrices. 
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
%   BEAM_PLANE_STD_DEV = Standard deviation (in pixels) of the intensity
%       distribution of the Gaussian-shaped beam in the Z-direction (out
%       of plane direction).
%       
%   TRANSFORMS = 4 x 4 x N array of N homogeneous matrices that
%   specify the affine transformations relating the [x, y, z] 
%   positions of the particles in the first image to the [x', y', z'] 
%   positions of the particles in the N'th image.
%
%   DIMS = Dimensionality of the images; DIMS = 2 generates 2-D images,
%       and DIMS = 3 generates 3-D volumes.
%
%   PARALLEL_PROCESSING = Boolean flag specifying whether or not
%       to run the code using parallel processing (1) or serial
%       processing (0).
%
%   RUN_COMPILED = Boolean flag specifying whether or not to run the
%       compiled version of the particle image generation code.
%       The compiled executes about 100 times faster than the
%       non-compiled code.
%
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
if nargin < 11
    COMPILED = 0;
end;

% Count the number of images to create
num_images = size(TRANSFORMS, 3);

% Pick between 2-D and 3-D images.
switch DIMS
    % Allocate the arrays for 2-D images
    case 2
        % Allocate the image list
        img_series_2D = zeros(IMAGE_HEIGHT, IMAGE_WIDTH,...
            num_images, 'uint16');
        
        % Just define img_series_3D so that the parfor loop runs
        img_series_3D = [];
        
    % Allocate the arrays for 3-D images.
    otherwise
        
        % Just define img_series_2D so that the parfor loop runs
        img_series_2D = [];
        
        % Allocate the image list
        img_series_3D = zeros(IMAGE_HEIGHT, IMAGE_WIDTH, ...
            IMAGE_DEPTH, num_images, 'uint16');

end

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

% Choose between serial and parallel processing
if PARALLEL_PROCESSING
    parfor k = 1 : num_images

            % Inform the user
            fprintf(1, 'Generating image %d of %d\n', k, num_images);

            % Transform the particle coordinates
            [Y2, X2, Z2] = transformImageCoordinates_3D(TRANSFORMS(:, :, k), ...
                X1, Y1, Z1, [yc, xc, zc]);

            % Set the particle max intensities to be proportional to the 
            % Gaussian intensity profile of the beam
            particleMaxIntensities = exp(-(Z2 - zc).^2 ./ ...
            (2 * BEAM_PLANE_STD_DEV ^ 2));

            % Pick between 2-D and 3-D images.
            switch DIMS
                case 2

                    % Generate the subsequent images.
                    if COMPILED
                        % Generate the second image (compiled code)
                        image_02 = generateParticleImage_mex( ...
                            augmentedHeight, augmentedWidth, X2, Y2, ...
                            particle_diameters, particleMaxIntensities);
                    else
                        % Generate the second image (non-compiled code)
                        image_02 = generateParticleImage_3D( ...
                            augmentedHeight, augmentedWidth, X2, Y2, ...
                            particle_diameters, particleMaxIntensities);
                    end

                    % Crop the second image
                    Image2Cropped =  flipud(image_02(...
                        augmentedHeight / 4 + 1 : ...
                        3 * augmentedHeight / 4, ...
                        augmentedWidth / 4 + 1 : ...
                        3 * augmentedWidth / 4));

                    % Recast as 16-bit
                    img_series_2D(:, :, k) = uint16( (2^16 - 1) .*...
                        Image2Cropped .*  2.8^2 ./ PARTICLE_DIAMETER_MEAN^2);

                otherwise

                    % Generate the second image
                    if COMPILED
                        % Run compiled MEX code
                        image_02 = generateParticleImage_3D_mex(...
                            augmentedHeight, augmentedWidth, ...
                            augmentedDepth, X2, Y2, Z2, ...
                            particle_diameters, particleMaxIntensities);
                    else
                        % Generate the second image
                        image_02 = generateParticleImage_3D(...
                            augmentedHeight, augmentedWidth, ...
                            augmentedDepth, X2, Y2, Z2, ...
                            particle_diameters, particleMaxIntensities);
                    end

                    % Crop the second image
                    Image2Cropped =  flipud(image_02(...
                        augmentedHeight / 4 + 1 : ...
                        3 * augmentedHeight / 4, ...
                        augmentedWidth / 4 + 1 : ...
                        3 * augmentedWidth / 4, ...
                        augmentedDepth / 4 + 1 : ...
                        3 * augmentedDepth / 4));


                % Recast as 16-bit
                img_series_3D(:, :, :, k) = uint16( (2^16 - 1) .* ...
                    Image2Cropped .*  2.8^2 ./ PARTICLE_DIAMETER_MEAN^2);

            end      
    end

else
    
    for k = 1 : num_images

            % Inform the user
            fprintf(1, 'Generating image %d of %d\n', k, num_images);

            % Transform the particle coordinates
            [Y2, X2, Z2] = transformImageCoordinates_3D(TRANSFORMS(:, :, k), ...
                X1, Y1, Z1, [yc, xc, zc]);

            % Set the particle max intensities to be proportional to the 
            % Gaussian intensity profile of the beam
            particleMaxIntensities = exp(-(Z2 - zc).^2 ./ ...
            (2 * BEAM_PLANE_STD_DEV ^ 2));

            % Pick between 2-D and 3-D images.
            switch DIMS
                case 2

                    % Generate the subsequent images.
                    if COMPILED
                        % Generate the second image (compiled code)
                        image_02 = generateParticleImage_mex( ...
                            augmentedHeight, augmentedWidth, X2, Y2, ...
                            particle_diameters, particleMaxIntensities);
                    else
                        % Generate the second image (non-compiled code)
                        image_02 = generateParticleImage_3D( ...
                            augmentedHeight, augmentedWidth, X2, Y2, ...
                            particle_diameters, particleMaxIntensities);
                    end

                    % Crop the second image
                    Image2Cropped =  flipud(image_02(...
                        augmentedHeight / 4 + 1 : ...
                        3 * augmentedHeight / 4, ...
                        augmentedWidth / 4 + 1 : ...
                        3 * augmentedWidth / 4));

                    % Recast as 16-bit
                    img_series_2D(:, :, k) = uint16( (2^16 - 1) .*...
                        Image2Cropped .*  2.8^2 ./ PARTICLE_DIAMETER_MEAN^2);

                otherwise

                    % Generate the second image
                    if COMPILED
                        % Run compiled MEX code
                        image_02 = generateParticleImage_3D_mex(...
                            augmentedHeight, augmentedWidth, ...
                            augmentedDepth, X2, Y2, Z2, ...
                            particle_diameters, particleMaxIntensities);
                    else
                        % Generate the second image
                        image_02 = generateParticleImage_3D(...
                            augmentedHeight, augmentedWidth, ...
                            augmentedDepth, X2, Y2, Z2, ...
                            particle_diameters, particleMaxIntensities);
                    end

                    % Crop the second image
                    Image2Cropped =  flipud(image_02(...
                        augmentedHeight / 4 + 1 : ...
                        3 * augmentedHeight / 4, ...
                        augmentedWidth / 4 + 1 : ...
                        3 * augmentedWidth / 4, ...
                        augmentedDepth / 4 + 1 : ...
                        3 * augmentedDepth / 4));

                    % Recast as 16-bit
                    % Matlab is giving a warning in this code
                    % that says the array img_series_3D is growing
                    % inside the loop, but in reality it will not grow
                    % iteratively because this line is only accessed
                    % when DIMS = 3, in which case img_series_3D will have
                    % been properly allocated.
                    img_series_3D(:, :, :, k) = uint16( (2^16 - 1) .* ...
                        Image2Cropped .*  2.8^2 ./ PARTICLE_DIAMETER_MEAN^2);

            end        
    end

end

% Assign the appropriate version of the image series
% to the output variable.
switch DIMS
    case 2
        IMAGE_SERIES = img_series_2D;
    otherwise
        IMAGE_SERIES = img_series_3D;
end

   


end % End of function

