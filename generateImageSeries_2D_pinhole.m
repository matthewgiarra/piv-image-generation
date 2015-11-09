function [imageMatrix1, imageMatrix2] = generateImageSeries_2D_pinhole(WORLD_DOMAIN, PARTICLE_DIAMETER_MEAN, ...
    PARTICLE_DIAMETER_STD, PARTICLE_CONCENTRATION, BEAM_PLANE_CENTER_WORLD_Z, BEAM_PLANE_STD_DEV, ...
    PARTICLE_TRANSFORMATION_MATRIX, CAMERA_PARAMETERS, RUN_COMPILED)
%[IMAGE_01, IMAGE_02] = generateImagePair_ND_mc(IMAGE_HEIGHT, ...
%     IMAGE_WIDTH, IMAGE_DEPTH, PARTICLE_DIAMETER_MEAN, ...
%     PARTICLE_DIAMETER_STD, PARTICLE_CONCENTRATION, BEAM_PLANE_STD_DEV, ...
%     PARTICLE_TRANSFORMATION_MATRIX, DIMS, RUN_COMPILED)
% Generates a series of 2-D or 3-D synthetic particle images
%   that have been transformed by a specified series of
%   transformation matrices. 
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
%   BEAM_PLANE_STD_DEV = Standard deviation (in pixels) of the intensity
%       distribution of the Gaussian-shaped beam in the Z-direction (out
%       of plane direction).
%       
%   PARTICLE_TRANSFORMATION_MATRIX = 4 x 4 homogeneous matrix that specifies the
%   affine transformation relating the [x, y, z] positions of the particles
%   in the first image to the [x', y', z'] particle positions in the
%   second image.
%
%   DIMS = Dimensionality of the images; DIMS = 2 generates 2-D images,
%       and DIMS = 3 generates 3-D volumes.
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

% Particle transforms
tforms = PARTICLE_TRANSFORMATION_MATRIX;

% Number of images to generate
num_images = size(tforms, 3);

% Count the number of cameras
num_cameras = CAMERA_PARAMETERS.NumberOfCameras;

% Image width and height in pixels
% Currently these are forced to all be the same so that they can be added
image_height_pixels = CAMERA_PARAMETERS.Cameras(1).Intrinsic.Pixel.Number.Rows;
image_width_pixels  = CAMERA_PARAMETERS.Cameras(1).Intrinsic.Pixel.Number.Columns;

% Double the height, width, and depth of the images so that rotations don't cause
% cropping.
world_domain_x = WORLD_DOMAIN(1, :);
world_domain_y = WORLD_DOMAIN(2, :);
world_domain_z = WORLD_DOMAIN(3, :);

% Domain limits
% X coordinate limits (l = lower limit, u = upper limit)
xl = world_domain_x(1);
xu = world_domain_x(2);
% Y coordinate limits (l = lower limit, u = upper limit)
yl = world_domain_y(1);
yu = world_domain_y(2);
% Z coordinate limits (l = lower limit, u = upper limit)
zl = world_domain_z(1);
zu = world_domain_z(2);

% Domain width, height, depth
domain_width  = xu - xl;
domain_height = yu - yl;
domain_depth  = zu - zl;

% Number of particles to generate
% Particle concentration in particles per cubic meter.
num_particles = round(PARTICLE_CONCENTRATION * ...
domain_width * domain_height * domain_depth);

% num_particles = 1000;

% Reference plane distance
L = CAMERA_PARAMETERS.TargetPlaneDistance;

% Randomly generate the world coordinates of particles
% Generate horizontal locations of particles in first image (column vector)
X1 = xl + (xu  - xl) * rand(num_particles, 1);
Y1 = yl + (yu  - yl) * rand(num_particles, 1);
Z1 = zl + (zu  - zl) * rand(num_particles, 1);

% Create a normal distribution of particle diameters
% Need to change this to be something physical!!!
particle_diameters = repmat(abs(PARTICLE_DIAMETER_STD * randn(num_particles, 1) ...
    + PARTICLE_DIAMETER_MEAN), [num_cameras, 1]);

	% Vectors of positions for debugging stuff
x_cam_01 = zeros(num_particles, num_cameras);
y_cam_01 = zeros(num_particles, num_cameras);
x_cam_02 = zeros(num_particles, num_cameras);
y_cam_02 = zeros(num_particles, num_cameras);

% Allocate the camera images
imageMatrix1 = zeros(image_height_pixels, image_width_pixels, 1);
imageMatrix2 = zeros(image_height_pixels, image_width_pixels, num_images);
		
% Set the particle max intensities to be proportional to the
% Gaussian intensity profile of the beam
particleMaxIntensities = repmat(exp(-(Z1 - BEAM_PLANE_CENTER_WORLD_Z).^2 ./ ...
    (2 * BEAM_PLANE_STD_DEV ^ 2)), [num_cameras, 1]);
	
% Generate the particle positions in the coordinate system of each camera
for c = 1 : num_cameras

	% Read the camera matrix
	camera_matrix = CAMERA_PARAMETERS.Cameras(c).Camera_Matrix;

	% Calculate image coordinates of each particle (pixel coordinates)
	[x_cam(:, c), y_cam(:, c)] = pinhole_camera_coordinate_transform(X1, Y1, Z1, camera_matrix);

end

% Render the images
if RUN_COMPILED

	% Generate the first image. Add the intensities to any already generated images.
	imageMatrix1 = generateParticleImage_mex(image_height_pixels, image_width_pixels, ...
		x_cam(:), y_cam(:), particle_diameters, particleMaxIntensities);

else

	% Generate the first image. Add the intensities to any already generated images.
	imageMatrix1 = generateParticleImage(image_height_pixels, image_width_pixels, ...
		x_cam(:), y_cam(:), particle_diameters, particleMaxIntensities);
	
end

% Image matrix 2
for n = 1 : num_images
	
	% Inform the user
	fprintf('Generating image %d of %d...\n', n, num_images);
	
	% Transform the particle coordinates
	[Y2, X2, Z2] = transformImageCoordinates_3D(tforms(:, :, n), ...
	    X1, Y1, Z1, [0, 0, 0]);
		
	% Set the particle max intensities to be proportional to the
	% Gaussian intensity profile of the beam
	particleMaxIntensities = repmat(exp(-(Z2 - BEAM_PLANE_CENTER_WORLD_Z).^2 ./ ...
	    (2 * BEAM_PLANE_STD_DEV ^ 2)), [num_cameras, 1]);
		
    % Generate the particle positions in the coordinate system of each camera
    for c = 1 : num_cameras

        % Read the camera matrix
        camera_matrix = CAMERA_PARAMETERS.Cameras(c).Camera_Matrix;

        % Calculate image coordinates of each particle (pixel coordinates)
        [x_cam(:, c), y_cam(:, c)] = pinhole_camera_coordinate_transform(X2, Y2, Z2, camera_matrix);

    end    
    
	
    % Render the images
    if RUN_COMPILED

        % Generate the first image. Add the intensities to any already generated images.
        imageMatrix2(:, :, n) = generateParticleImage_mex(image_height_pixels, image_width_pixels, ...
        x_cam(:), y_cam(:), particle_diameters, particleMaxIntensities);

    else

        % Generate the first image. Add the intensities to any already generated images.
        imageMatrix2(:, :, n) = generateParticleImage(image_height_pixels, image_width_pixels, ...
        x_cam(:), y_cam(:), particle_diameters, particleMaxIntensities);

    end
   
	
end

% Make a plot of particle positions
% plot_colors = 'krbgcmy';
% Plot marker
% plot_marker = '.';
% figure(1);
% subplot(1, 2, 1);
% hold off;
% for k = 1 : num_cameras
%
% 	% Plot marker string
% 	plot_marker_string = [plot_marker, plot_colors(k)];
%
% 	% Plot the particle positions
% 	plot(x_vect_01(:, k), y_vect_01(:, k), plot_marker_string);
%
% 	% Hold plot
% 	hold on;
%
% end
%
% hold off;
% axis image;
%
% subplot(1, 2, 2);
% hold off;
% for k = 1 : num_cameras
%
% 	% Plot marker string
% 	plot_marker_string = [plot_marker, plot_colors(k)];
%
% 	% Plot the particle positions
% 	plot(x_vect_02(:, k), y_vect_02(:, k), plot_marker_string);
%
% 	% Hold plot
% 	hold on;
%
% end
% hold off
% axis image

% Convert the first image to 16 bit and save it to the output variable.
imageMatrix1 = ( (2^16 - 1) * imageMatrix1 * ...
2.8 ^ 2 / PARTICLE_DIAMETER_MEAN ^ 2);

% Convert the second image to 16 bit and save it to the output variable.
imageMatrix2 = ( (2^16 - 1) * imageMatrix2 * ...
    2.8 ^ 2 / PARTICLE_DIAMETER_MEAN ^ 2);


end % End of function

