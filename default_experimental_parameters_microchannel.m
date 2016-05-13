function OUTPUT_PARAMS = default_experimental_parameters_microchannel(INPUT_PARAMS);

% Channel size
channel_width_microns = 5E3;
channel_depth_microns = 100;

% Target displacement in pixels
dx_target_pix = 15.0;

% Temperature
T_kelvin = 300;

% Viscosity in Pa * S
viscosity_pas = 1.12E-3;

% Flow rate in µl/min
flow_rate_ul_min = 5;

% Objective magnification
objective_magnification = 60;

% Particle diameter in microns
dp_microns = 0.1;

% Stock particle solution volume fraction
% Particles come as 1% solids
volume_fraction_particles_stock = 1E-2;

% Particle dilution factor
% This is the factor by which 
% we dilute the particle suspension
particle_dilution_factor = 25;

% Volume fraction of particles
volume_fraction_particles = volume_fraction_particles_stock / ...
    particle_dilution_factor;

% Pixel size in microns
pixel_size_um = INPUT_PARAMS.Image.PixelSize;

% Temperature in kelvin
OUTPUT_PARAMS.Temperature = T_kelvin;

% Viscosity in Pa * S
OUTPUT_PARAMS.Viscosity = viscosity_pas;

% Channel dimensions in microns
OUTPUT_PARAMS.ChannelWidth = channel_width_microns;
OUTPUT_PARAMS.ChannelDepth = channel_depth_microns; 

% Flow rate in µl/minute
% I chose these units to be 
% consistent with the syringe pump.
OUTPUT_PARAMS.FlowRate = 5;

% Objective magnification 
OUTPUT_PARAMS.Objective = objective_magnification;

% Range of particle concentrations (Volume fraction particles)
OUTPUT_PARAMS.ParticleConcentration = volume_fraction_particles * [1, 1];

% Particle diameter in microns
OUTPUT_PARAMS.ParticleDiameter = dp_microns;

% Diffusion
diffusion_std_dev_pix = calculate_particle_diffusion(...
    channel_width_microns, channel_depth_microns, ...
    flow_rate_ul_min, objective_magnification, ...
    pixel_size_um, T_kelvin, dp_microns, viscosity_pas, dx_target_pix); 

% Diffusion
OUTPUT_PARAMS.DiffusionStdDev =  diffusion_std_dev_pix * [1, 1];

end