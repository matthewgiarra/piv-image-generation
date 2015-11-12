function IMAGE = generate_micro_piv_image(X_PIX, Y_PIX, Z_PIX, Parameters)

% Parse structure: image dimensions in pixels
region_width_pixels  = Parameters.Image.Width;
region_height_pixels = Parameters.Image.Height;

% Image magnification (microns per pixel)
pixel_size_microns = Parameters.Image.PixelSize;

% Optics parameters
%
% Objective lens
% Magnification (unitless)
objective_magnification = Parameters.Optics.Objective.Magnification;

% Objective lens focal length in microns
focal_length_microns = Parameters.Optics.Objective.FocalLength;

% Objective lens numerical aperture (unitless)
NA = Parameters.Optics.Objective.NA;

% Laser wavelength in microns
wavelength_microns = Parameters.Optics.Laser.Wavelength;

% Fraction above background intensity to render particles
intensity_fraction = Parameters.Experiment.IntensityFraction;

% Particle diameter
% Particle diameter in microns
particle_diameter_microns = Parameters.Experiment.ParticleDiameter;

% Channel depth in microns
channel_depth_microns = Parameters.Experiment.ChannelDepth;

% Particle concentration (particles per µm^3)
particle_concentration = Parameters.Experiment.ParticleConcentration;

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

% Render the first image.
IMAGE = generateParticleImage(...
    region_height_pixels, region_width_pixels, ...
    X_render, Y_render, ...
    particle_image_diameters_pixels_render, ...
    particle_max_intensities_render) + background_intensity;

end