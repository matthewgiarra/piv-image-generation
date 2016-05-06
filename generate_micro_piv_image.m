function IMAGE = generate_micro_piv_image(X_PIX, Y_PIX, Z_PIX, ...
    region_height_pixels, region_width_pixels, particle_diameter_microns ,...
    pixel_size_microns, particle_concentration, channel_depth_microns, ...
    objective_magnification, NA, working_distance_microns, focal_length_microns,...
    wavelength_microns, intensity_fraction)
    
    
% Limit of the near-field in microns
% This is calculated from the description 
% in Eckstein & Vlachos 2009 "Digital particle image velocimetry (DPIV)
% robust phase correlation"
% by setting the size of a particle (in microns) at 
% a distance z from the focal plane (taken from equation 21) 
% equal to twice the average particle spacing (concentration in particles /
% / microns^3) ^(-1/3)
% This is here is half the depth of field -- the absolute distance from
% the focal plane corresponding to the in field depths.
depth_of_field_microns = sqrt(4 * particle_concentration^(-2/3) - particle_diameter_microns^2 - 5.95 * wavelength_microns^2 / (4 * NA^2));

% Number of particles
num_particles = length(X_PIX(:));

% Vector of particle diameters
dp = particle_diameter_microns * ones(num_particles, 1);

% Convert Z position to microns
Z_um = Z_PIX * pixel_size_microns / objective_magnification;

% % Temporary
focal_plane_z_microns = channel_depth_microns /2 ;

% Z position in the coordinate system of the focal plane
z_um_focal = Z_um - focal_plane_z_microns;

% Particles that lie within the depth of field
inds = find(abs(z_um_focal) < depth_of_field_microns);

% Calculate particle max intensities and diameters
% These are in units of microns
[particle_max_intensities, particle_image_diameters_microns]  = ...
    calculate_particle_max_intensity(objective_magnification, ...
    dp(inds), wavelength_microns, NA, z_um_focal(inds));

% Convert particle image diameters to pixels
particle_image_diameters_pixels = ...
    particle_image_diameters_microns / pixel_size_microns;

% Object plane
object_plane_microns = channel_depth_microns / 2;

% Background intensity
background_intensity = calculate_background_intensity(...
    channel_depth_microns, particle_concentration, ...
    particle_diameter_microns, wavelength_microns,...
    objective_magnification, NA, working_distance_microns, ...
    focal_length_microns, object_plane_microns);

% Render the image
IMAGE = generateParticleImage(...
    region_height_pixels, region_width_pixels, ...
    X_PIX(inds), Y_PIX(inds), ...
    particle_image_diameters_pixels, ...
    particle_max_intensities) + background_intensity;

end

