function  [IMAGE1, IMAGE2] = generateImagePair_micro_mc(Parameters)

% % % % % % % %  % 
% BEGIN FUNCTION %
% % % % % % % %  % 

% Parse structure: image dimensions in pixels
region_width_pixels  = Parameters.Image.Width;
region_height_pixels = Parameters.Image.Height;

% Image magnification (microns per pixel)
pixel_size_microns = Parameters.Image.PixelSize;

% Optics parameters
objective_magnification = Parameters.Optics.Objective.Magnification;
%
% Image transformation
transformation_matrix = Parameters.Transformation;

% Experiment parameters
% 
% Channel depth in microns
channel_depth_microns = Parameters.Experiment.ChannelDepth;

% Particle concentration (particles per µm^3)
particle_concentration = Parameters.Experiment.ParticleConcentration;

% Diffusion
particle_diffusion_stdev = Parameters.Experiment.DiffusionStDev;

% Size of the X domain in microns
domain_x_um = (2 * region_width_pixels + 1) * ...
    pixel_size_microns / objective_magnification;

% Size of the Y domain in microns
domain_y_um = (2 * region_height_pixels + 1) * ...
    pixel_size_microns / objective_magnification;

% Size of the Z domain in microns
domain_z_um = channel_depth_microns;

% Size of the Z domain in pixels
region_depth_pixels = domain_z_um * ...
    objective_magnification / pixel_size_microns;

% Number of particles to generate
nParticles = round(particle_concentration ...
    * domain_x_um * domain_y_um * domain_z_um);

% Randomly generate the coordinates of particles
% Generate horizontal locations of particles in first image (column vector)
X1_pix = region_width_pixels  * (2 * rand(nParticles, 1) - 1);
Y1_pix = region_height_pixels * (2 * rand(nParticles, 1) - 1);
Z1_pix = region_depth_pixels  * (2 * rand(nParticles, 1) - 1);

% Generate the first image
IMAGE1 = generate_micro_piv_image(X1_pix, Y1_pix, Z1_pix, Parameters);

% Transform particle coordinates
[X2_pix, Y2_pix, Z2_pix] = transformImageCoordinates_3D(...
    X1_pix, Y1_pix, Z1_pix, transformation_matrix);

% Add diffusion
X2_pix = X2_pix + particle_diffusion_stdev * randn(nParticles, 1);
Y2_pix = Y2_pix + particle_diffusion_stdev * randn(nParticles, 1);
Z2_pix = Z2_pix + particle_diffusion_stdev * randn(nParticles, 1);

% Generate the first image
IMAGE2 = generate_micro_piv_image(X2_pix, Y2_pix, Z2_pix, Parameters);


end % End of function

