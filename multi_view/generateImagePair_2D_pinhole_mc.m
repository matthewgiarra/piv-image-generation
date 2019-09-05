function [IMAGE_01, IMAGE_02] = generateImagePair_2D_pinhole_mc(WORLD_DOMAIN, PARTICLE_DIAMETER_MEAN, ...
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

% Double the height, width, and depth of the images so that rotations don't cause
% cropping.
world_domain_x = WORLD_DOMAIN(1, :);
world_domain_y = WORLD_DOMAIN(2, :);
world_domain_z = WORLD_DOMAIN(3, :);

% Count the number of cameras
num_cameras = CAMERA_PARAMETERS.NumberOfCameras;

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

% Randomly generate the world coordinates of particles
% Generate horizontal locations of particles in first image (column vector)
X1 = xl + (xu  - xl) * rand(num_particles, 1);
Y1 = yl + (yu  - yl) * rand(num_particles, 1);
Z1 = zl + (zu  - zl) * rand(num_particles, 1);

% Transform the particle coordinates
[Y2, X2, Z2] = transformImageCoordinates_3D(PARTICLE_TRANSFORMATION_MATRIX, ...
    X1, Y1, Z1, [0, 0, 0]);
	
% Create a normal distribution of particle diameters
% Need to change this to be something physical!!!
particle_diameters = repmat(abs(PARTICLE_DIAMETER_STD * randn(num_particles, 1) ...
    + PARTICLE_DIAMETER_MEAN), [num_cameras, 1]);

% Gaussian Function that expresses the intensity distribution of the 
% particle image on the "sensor" This is set to 1 for now since the 
% depth-positions of the particles are known for 3D rendering.
% particleMaxIntensities = ones(num_particles, 1);

% % Set the particle max intensities to be proportional to the
% % Gaussian intensity profile of the beam
particleMaxIntensities_01 = repmat(exp(-(Z1 - BEAM_PLANE_CENTER_WORLD_Z).^2 ./ ...
    (2 * BEAM_PLANE_STD_DEV ^ 2)), [num_cameras, 1]);

% Set the particle max intensities to be proportional to the
% Gaussian intensity profile of the beam
particleMaxIntensities_02 = repmat(exp(-(Z2 - BEAM_PLANE_CENTER_WORLD_Z).^2 ./ ...
    (2 * BEAM_PLANE_STD_DEV ^ 2)), [num_cameras, 1]);


% Fix the brightness
% particleMaxIntensities_01 = ones(num_cameras * num_particles, 1);
% particleMaxIntensities_02 = ones(num_cameras * num_particles, 1);
	
% Image width and height in pixels
% Currently these are forced to all be the same so that they can be added
image_height_pixels = CAMERA_PARAMETERS.Cameras(1).Intrinsic.Pixel.Number.Rows;
image_width_pixels  = CAMERA_PARAMETERS.Cameras(1).Intrinsic.Pixel.Number.Columns;

% Vectors of positions for debugging stuff
x_vect_01 = zeros(num_particles, num_cameras);
y_vect_01 = zeros(num_particles, num_cameras);
x_vect_02 = zeros(num_particles, num_cameras);
y_vect_02 = zeros(num_particles, num_cameras);

% Make a plot of particle positions
plot_colors = 'krbgcmy';

% Plot marker
plot_marker = '.';


% Generate the particle positions in the coordinate system of each camera
for k = 1 : num_cameras
	
	% Read the camera matrix
	camera_matrix = CAMERA_PARAMETERS.Cameras(k).Camera_Matrix;
	
	% Calculate image coordinates of each particle (pixel coordinates)
	[x_cam_01, y_cam_01] = pinhole_camera_coordinate_transform(X1, Y1, Z1, camera_matrix);
	[x_cam_02, y_cam_02] = pinhole_camera_coordinate_transform(X2, Y2, Z2, camera_matrix);
	
	x_vect_01(:, k) = x_cam_01;
	y_vect_01(:, k) = y_cam_01;
	
	x_vect_02(:, k) = x_cam_02;
	y_vect_02(:, k) = y_cam_02;
	
end

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
	
% Render the images
if RUN_COMPILED
	
	% Generate the first image. Add the intensities to any already generated images.
	particle_image_01 = generateParticleImage_mex(image_height_pixels, image_width_pixels, ...
		x_vect_01(:), y_vect_01(:), particle_diameters, particleMaxIntensities_01);
	
	% Generate the first image. Add the intensities to any already generated images.
	particle_image_02 = generateParticleImage_mex(image_height_pixels, image_width_pixels, ...
		x_vect_02(:), y_vect_02(:), particle_diameters, particleMaxIntensities_02);
	
else
	
	% Generate the first image. Add the intensities to any already generated images.
	particle_image_01 = generateParticleImage(image_height_pixels, image_width_pixels, ...
		x_vect_01(:), y_vect_01(:), particle_diameters, particleMaxIntensities_01);
	
	% Generate the first image. Add the intensities to any already generated images.
	particle_image_02 = generateParticleImage(image_height_pixels, image_width_pixels, ...
		x_vect_02(:), y_vect_02(:), particle_diameters, particleMaxIntensities_02);
		
end
	


% Convert the first image to 16 bit and save it to the output variable.
IMAGE_01 = ( (2^16 - 1) .* particle_image_01 .* ...
2.8 ^ 2 / PARTICLE_DIAMETER_MEAN ^ 2);

% Convert the second image to 16 bit and save it to the output variable.
IMAGE_02 = ( (2^16 - 1) .* particle_image_02 .* ...
    2.8 ^ 2 ./ PARTICLE_DIAMETER_MEAN ^ 2);

end % End of function

