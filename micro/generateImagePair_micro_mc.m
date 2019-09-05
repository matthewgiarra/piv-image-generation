function  [IMAGE1, IMAGE2] = generateImagePair_micro_mc(...
    region_height_pixels, region_width_pixels, particle_diameter_microns ,...
    pixel_size_microns, particle_concentration, channel_depth_microns, ...
    objective_magnification, NA, working_distance_microns, focal_length_microns,...
    wavelength_microns, diffusion_std_dev, Tform)

% % % % % % % %  % 
% BEGIN FUNCTION %
% % % % % % % %  % 

% Size of the X domain in microns
domain_x_um = (2 * region_width_pixels + 1) * ...
    pixel_size_microns / objective_magnification;

% Size of the Y domain in microns
domain_y_um = (2 * region_height_pixels + 1) * ...
    pixel_size_microns / objective_magnification;

% Size of the near-field Z domain in microns
domain_z_um = channel_depth_microns;

% Size of the Z domain in pixels
region_depth_pixels = domain_z_um * ...
    objective_magnification / pixel_size_microns;

% Volume of domain in µm^3
domain_volume = domain_x_um * domain_y_um * domain_z_um;

% Number of particles to generate
% nParticles = round(particle_concentration ...
%     * domain_x_um * domain_y_um * domain_z_um);

nParticles = round(6 * particle_concentration * domain_volume ...
    / (pi * particle_diameter_microns^3));

% Randomly generate the coordinates of particles
% Generate horizontal locations of particles in first image (column vector)
X1_pix = region_width_pixels  * (2 * rand(nParticles, 1) - 1);
Y1_pix = region_height_pixels * (2 * rand(nParticles, 1) - 1);
Z1_pix = region_depth_pixels  * rand(nParticles, 1);

% Transform particle coordinates
[X2_pix, Y2_pix, Z2_pix] = transformImageCoordinates_3D(...
    X1_pix, Y1_pix, Z1_pix, Tform);

% Add diffusion
X2_pix = X2_pix + diffusion_std_dev * randn(nParticles, 1);
Y2_pix = Y2_pix + diffusion_std_dev * randn(nParticles, 1);
Z2_pix = Z2_pix + diffusion_std_dev * randn(nParticles, 1);

% u = Tform(1, 4);
% v = Tform(2, 4);
% Skip = 100;
% x1 = X1_pix(1 : Skip: end);
% y1 = Y1_pix(1 : Skip: end);
% x2 = X2_pix(1 : Skip: end);
% y2 = Y2_pix(1 : Skip: end);
% dx = x2 - x1;
% dy = y2 - y1;
% plot(x1, y1, 'ok', 'markerfacecolor', 'black');
% hold on
% quiver(x1, y1, dx, dy, 0, 'black');
% hold off
% axis image;
% ylim([1, 127]);
% xlim([1, 127]);
% title({sprintf('u = %0.2f, v = %0.2f pix', u, v), ...
%     sprintf('Diffusion std dev: %0.2f pix', diffusion_std_dev)}...
%     , 'FontSize', 20);
% set(gca, 'FontSize', 12)
% box on;

% Generate the first image
IMAGE1 = generate_micro_piv_image(X1_pix, Y1_pix, Z1_pix, ...
    region_height_pixels, region_width_pixels, particle_diameter_microns ,...
    pixel_size_microns, particle_concentration, channel_depth_microns, ...
    objective_magnification, NA, working_distance_microns, focal_length_microns,  ...
    wavelength_microns);


% Generate the second image
% Inform the user
IMAGE2 = generate_micro_piv_image(X2_pix, Y2_pix, Z2_pix, ...
    region_height_pixels, region_width_pixels, particle_diameter_microns ,...
    pixel_size_microns, particle_concentration, channel_depth_microns, ...
    objective_magnification, NA, working_distance_microns, focal_length_microns, ...
    wavelength_microns);

% Max value of the two images
% Both images will be scaled
% according to this value
% so that their exposures
% are the same.
max_val = max([IMAGE1(:); IMAGE2(:)]);

% Scale the images by the max value.
IMAGE1 = 0.95 * IMAGE1 ./ max_val;
IMAGE2 = 0.95 * IMAGE2 ./ max_val;

end % End of function










