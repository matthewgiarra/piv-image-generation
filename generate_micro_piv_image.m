function IMAGE = generate_micro_piv_image(X_PIX, Y_PIX, Z_PIX, ...
    region_height_pixels, region_width_pixels, particle_diameter_microns ,...
    pixel_size_microns, particle_concentration, channel_depth_microns, ...
    objective_magnification, focal_length_microns, NA, ...
    wavelength_microns, intensity_fraction)
    
    
% Number of particles
num_particles = length(X_PIX(:));

% Vector of particle diameters
dp = particle_diameter_microns * ones(num_particles, 1);

% Convert Z position to microns
Z_um = Z_PIX * pixel_size_microns / objective_magnification;

% Calculate particle max intensities and diameters
% These are in units of microns
[particle_max_intensities, particle_image_diameters_microns, ...
    emission_power]  = ...
    calculate_particle_max_intensity(objective_magnification, ...
    dp, wavelength_microns, NA, ...
    focal_length_microns, Z_um);

% Convert particle image diameters to pixels
particle_image_diameters_pixels = ...
    particle_image_diameters_microns / pixel_size_microns;

% Background intensity
background_intensity = calculate_background_intensity(...
    channel_depth_microns, particle_concentration,...
    objective_magnification, NA, emission_power);

% Indices of particles to render
render_inds = find(particle_max_intensities > ...
    intensity_fraction * background_intensity);

% Select the diameters of the particles to be rendered.
particle_image_diameters_pixels_render = ...
    particle_image_diameters_pixels(render_inds);

% Select the intensities of the particles to be rendered.
particle_max_intensities_render = particle_max_intensities(render_inds);

% Center pixels
xc = region_width_pixels  / 2 + 1;
yc = region_height_pixels / 2 + 1;

% Select the particles to be rendered.
X_render = X_PIX(render_inds) + xc;
Y_render = Y_PIX(render_inds) + yc;

% Render the image
IMAGE = generateParticleImage(...
    region_height_pixels, region_width_pixels, ...
    X_render, Y_render, ...
    particle_image_diameters_pixels_render, ...
    particle_max_intensities_render);

end

