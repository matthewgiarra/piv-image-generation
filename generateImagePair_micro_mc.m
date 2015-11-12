function  [IMAGE1, IMAGE2] = generateImagePair_micro_mc(...
    ...
    Parameters, RUN_COMPILED)


% % % % % % % % % 
% BEGIN FUNCTION %
% % % % % % % % % 

% Parse structure: image dimensions in pixels
region_width_pixels = Parameters.Image.Width;
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

% Image transformation
transformation_matrix = Parameters.Transformation;

% Experiment parameters
% 
% Channel depth in microns
channel_depth_microns = Parameters.Experiment.ChannelDepth;

% Particle diameter in microns
particle_diameter_microns = Parameters.Experiment.ParticleDiameter;

% Particle concentration (particles per µm^3)
particle_concentration = Parameters.Experiment.ParticleConcentration;

% Diffusion
particle_diffusion_stdev = Parameters.Experiment.DiffusionStDev;

% Domain limits
x_min_pix = -region_width_pixels;
x_max_pix =  region_width_pixels;
y_min_pix = -region_height_pixels;
y_max_pix =  region_height_pixels;

domain_x_um = (x_max_pix - x_min_pix) * pixel_size_microns / objective_magnification;
domain_y_um = (y_max_pix - y_min_pix) * pixel_size_microns / objective_magnification;
z_min_um = -channel_depth_microns / 2;
z_max_um =  channel_depth_microns / 2;
domain_z_um = z_max_um - z_min_um;

% Number of particles to generate
nParticles = round(particle_concentration ...
    * domain_x_um * domain_y_um * domain_z_um);

% Make random displacement matrix.
random_displacement_matrix = particle_diffusion_stdev * randn([nParticles, 2]); 

% Randomly generate the coordinates of particles
% Generate horizontal locations of particles in first image (column vector)
X1 = x_min_pix + (x_max_pix - x_min_pix) * rand(nParticles, 1);
Y1 = y_min_pix + (y_max_pix - y_min_pix) * rand(nParticles, 1);
Z1 = z_min_um + (z_max_um - z_min_um) * rand(nParticles, 1);
% Z1 = zeros(nParticles, 1);

% dp = abs(particle_diameter_microns + 1 * randn(nParticles, 1));

dp = particle_diameter_microns * ones(nParticles, 1);

% Calculate particle max intensities and diameters
% These are in units of microns
[particle_max_intensities, particle_image_diameters_microns, ...
    emission_power]  = ...
    calculate_particle_max_intensity(objective_magnification, ...
    dp, wavelength_microns, NA, ...
    focal_length_microns, Z1);

% emission_power = 1;

% emission_power_scaled = 0.3 * mean(emission_power);

% Background intensity
background_intensity = calculate_background_intensity(...
    channel_depth_microns, particle_concentration,...
    objective_magnification, NA, emission_power);

% Indices of particles to render
render_inds = find(particle_max_intensities > ...
    intensity_fraction * background_intensity);

% Convert particle image diameters to pixels
particle_image_diameters_pixels = ...
    particle_image_diameters_microns / pixel_size_microns;

% Particle image diameters in microns to render
particle_image_diameters_microns_render = particle_image_diameters_microns(render_inds);

% Coordinates of particles to render
X1_render = X1(render_inds) + region_width_pixels + 1;
Y1_render = Y1(render_inds) + region_height_pixels + 1;
Z1_render = Z1(render_inds);
particle_image_diameters_pixels_render = ...
    particle_image_diameters_pixels(render_inds);
particle_max_intensities_render = particle_max_intensities(render_inds);


% hold off
% for k = 1 : length(X1_render)
% plot(X1_render(k), Y1_render(k), '.k', 'markersize', particle_image_diameters_pixels_render(k));
% hold on;
% end
% hold off;
tic
% Generate the first image. Choose between compiled and scripted code.
if RUN_COMPILED
    image1 = generateParticleImage_mex(region_height_pixels, region_width_pixels, X1_render, Y1_render, particle_image_diameters_pixels_render, particle_max_intensities_render);
else
    image1 =     generateParticleImage(region_height_pixels, region_width_pixels, X1_render, Y1_render, particle_image_diameters_pixels_render, particle_max_intensities_render);
end
toc

imagesc(image1 + background_intensity); axis image;
colormap gray
c = caxis;
caxis([0 1 * max(c)]);


% Crop the image and flip it vertically to place it in a right-handed coordinate system.
Image1Cropped = flipud(image1(augmentedHeight / 4 + 1 : 3 * augmentedHeight / 4, augmentedWidth / 4 + 1 : 3 * augmentedWidth / 4));

% Make 16 bit. Save to output variable.
IMAGE1 = uint16( (2^16 - 1) .* Image1Cropped .* 2.8^2 / PARTICLE_DIAMETER_MEAN ^2);

% Transform the particle coordinates
[Y2, X2] = transformImageCoordinates(TRANSFORMATION_MATRIX, X1, Y1, [yc xc]);

% Generate the second image. Choose between compiled and scripted code.
if RUN_COMPILED
    image2 = generateParticleImage_mex(augmentedHeight, augmentedWidth, X2 + random_displacement_matrix(:, 2), Y2 + random_displacement_matrix(:, 1), particle_image_diameters, particle_max_intensities);
else
    image2 = generateParticleImage(    augmentedHeight, augmentedWidth, X2 + random_displacement_matrix(:, 2), Y2 + random_displacement_matrix(:, 1), particle_image_diameters, particle_max_intensities);
end

% Crop the second image
Image2Cropped =  flipud(image2(augmentedHeight / 4 + 1: 3 * augmentedHeight / 4, augmentedWidth / 4 + 1: 3 * augmentedWidth / 4));

% Recast as 16-bit
IMAGE2 = uint16( (2^16 - 1) .* Image2Cropped .*  2.8^2 ./ PARTICLE_DIAMETER_MEAN^2);

end % End of function

