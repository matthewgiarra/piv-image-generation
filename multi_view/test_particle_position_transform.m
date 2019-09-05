function test_particle_position_transform(Camera_Parameters);

% Seed the number generator
rng(1);

% Number of particles
n_particles = 3E3;

% Particle diameter
particle_diameters = 2 * sqrt(8) * ones(n_particles, 1);

% Particle brightness
particle_max_intensity = 1 * ones(n_particles, 1);

% Number of pixels
n_pixels_rows = 128;
n_pixels_cols = 128;

% World axis limits
x_world_limits = 1 * [-1, 1];
y_world_limits = 1 * [-1, 1];
z_world_limits = 0.05 * [-1, 1];

% Count the number of cameras
% Long line because of input checking.
num_cameras = length(Camera_Parameters);

% X world limits (l = lower, u = upper)
xl = x_world_limits(1);
xu = x_world_limits(2);

% Y world limits (l = lower, u = upper)
yl = y_world_limits(1);
yu = y_world_limits(2);

% Z world limits (l = lower, u = upper)
zl = z_world_limits(1);
zu = z_world_limits(2);

% Discrete positions
x = xl + (xu - xl) * rand(n_particles, 1);
y = yl + (yu - yl) * rand(n_particles, 1);
z = zl + (zu - zl) * rand(n_particles, 1);

% Pixel size in meters (about 17 microns for Photrons)
pixel_size_x = 1.7E-5;
pixel_size_y = 1.7E-5;

% Plot colors
colors = 'rbkgcmy';

% Allocations
% Pixel sizes
pixel_size_x = zeros(num_cameras, 1);
pixel_size_y = zeros(num_cameras, 1);

% Sensor sizes
sensor_height = zeros(num_cameras, 1);
sensor_width  = zeros(num_cameras, 1);

% Calculate camera calibration matrices
camera_matrices = zeros(4, 4, num_cameras);
for k = 1 : num_cameras
	
	% Sensor size
	sensor_height(k) = Camera_Parameters(k).Intrinsic.Sensor.Height;
	sensor_width(k)  = Camera_Parameters(k).Intrinsic.Sensor.Width;
	
	% Camera rotation
	rx = Camera_Parameters(k).Extrinsic.Rotation.X;
	ry = Camera_Parameters(k).Extrinsic.Rotation.Y;
	rz = Camera_Parameters(k).Extrinsic.Rotation.Z;
	
	% Camera position (translation)
	tx = Camera_Parameters(k).Extrinsic.Translation.X;
	ty = Camera_Parameters(k).Extrinsic.Translation.Y;
	tz = Camera_Parameters(k).Extrinsic.Translation.Z;
	
	% Camera focal length
	f = Camera_Parameters(k).Intrinsic.FocalLength;
	
	% Pixel size
	pixel_size_x(k) = Camera_Parameters(k).Intrinsic.Pixel.Width;
	pixel_size_y(k) = Camera_Parameters(k).Intrinsic.Pixel.Height;
	
	% Pixel size vector
	pixel_size_vect = [pixel_size_y(k); pixel_size_x(k)];
	
	% Sensor size vector
	sensor_size_vect = [sensor_height(k); sensor_width(k)];
			
	% Calculate camera matrix (pixel coordinates)
	camera_matrices(:, :, k) = ...
		calculate_camera_matrix_pixel(rx, ry, rz, tx, ty, tz, f, sensor_size_vect, pixel_size_vect);
	
end

% Allocate vector particle positions for debugging
x_vect = zeros(n_particles, num_cameras);
y_vect = zeros(n_particles, num_cameras);

% Allocate the particle image array
particle_image = zeros(n_pixels_rows, n_pixels_cols);

figure(1);
% Loop over all the cameras
for k = 1 : num_cameras

	% Calculate image coordinates
	[x_cam, y_cam] = pinhole_camera_coordinate_transform(x, y, z, camera_matrices(:, :, k));

	% Plot color
	plot_color = colors(mod(k, length(colors)));

	% Plot color code
	color_code = ['.' plot_color];

	% Make the plot
	plot(x_cam, y_cam, color_code);
	
	% Vector positions for debugging stuff
	x_vect(:, k) = x_cam;
	y_vect(:, k) = y_cam;
	
	% Hold the plot
	hold on;
	
	% Add to the particle array
	particle_image = particle_image + ...
	 generateParticleImage(n_pixels_rows, n_pixels_cols, ...
	 x_cam(:), y_cam(:), ...
	 particle_diameters, particle_max_intensity);
end

% Format axes
axis image
set(gca, 'ydir', 'reverse');
xlim([0, sensor_width(1) / pixel_size_x(1)]);
ylim([0, sensor_height(1) / pixel_size_y(1)]);
% Release the plot hold
hold off;


figure(2);
imagesc(particle_image); axis image; colormap gray;



end




