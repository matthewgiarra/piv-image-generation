  function JOBLIST = MonteCarloImageGenerationJobFile_micro()

% Region dimensions
region_width_pixels  = 1024;
region_height_pixels  = 1024;

% Translations in pixels
% tx_pix = 7.5;
tx_pix = 15;
ty_pix = 0;

% Objective magnification
objective_magnification = 20;

% Pixel size in microns
pixel_size_microns = 10;

% Particle diameter in microns
dp_microns = 0.05;

% Particle concentration (volume fraction)
particle_volume_fraction = 1E-6;
 
DefaultJob.JobOptions.ParallelProcessing = 0;
DefaultJob.JobOptions.NumberOfDigits = 6;
DefaultJob.JobOptions.ReSeed = 1;

DefaultJob.ImageType = 'synthetic';
DefaultJob.SetType = 'mc';
DefaultJob.CaseName = 'piv_test_images';
% DefaultJob.ProjectRepository = '/home/shannon/b/aether/Projects/APC/';
DefaultJob.ProjectRepository = '~/Desktop/piv_test_images';

DefaultJob.Parameters.Image.Height = region_height_pixels;
DefaultJob.Parameters.Image.Width = region_width_pixels;
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1;
DefaultJob.Parameters.Sets.ImagesPerSet = 1;

% Rigid-body displacements (pixels)
DefaultJob.Parameters.Translation.X =  tx_pix * [1, 1];
DefaultJob.Parameters.Translation.Y =  ty_pix * [1, 1];
DefaultJob.Parameters.Translation.Z =  0  * [0, 2];

% Range of isotropic scaling factors
DefaultJob.Parameters.Scaling = 1 * [1 1]; 

% Range of rotation angles (degrees)
DefaultJob.Parameters.Rotation.Z_01 = 0 * [1 1];
DefaultJob.Parameters.Rotation.Y    = 0 * [1 1];
DefaultJob.Parameters.Rotation.Z_02 = 0 * [1 1];

% Shearing
DefaultJob.Parameters.Shear.XY = 0.0 * [1, 1];
DefaultJob.Parameters.Shear.XZ = 0.0 * [0, 1];
DefaultJob.Parameters.Shear.YX = 0 * [1, 1];
DefaultJob.Parameters.Shear.YZ = 0 * [0.00, 0.10];
DefaultJob.Parameters.Shear.ZX = 0 * [0.00, 0.10];
DefaultJob.Parameters.Shear.ZY = 0 * [0.00, 0.10];

% Noise parameters
DefaultJob.Parameters.Noise.Mean = 0.00 * [0, 0];
DefaultJob.Parameters.Noise.Std = 3E-2 * [1, 1];

% Optics
% Size of the pixels in microns
% Assumes square pixels
DefaultJob.Parameters.Image.PixelSize = pixel_size_microns;

% Bit depth of the imaging system
% Can be anything but typically 
% 8, 12, or 16
DefaultJob.Parameters.Image.BitDepth = 8;

% Objective lens parameters
DefaultJob.Parameters.Optics.Objective.Magnification =...
    objective_magnification;

% Laser wavelength in microns
DefaultJob.Parameters.Optics.Laser.Wavelength = 0.532;

% Flow rate in ul / min
DefaultJob.Parameters.Experiment.FlowRate = 5;

% Particle size in microns
DefaultJob.Parameters.Experiment.ParticleDiameter = dp_microns * [1, 1];
DefaultJob.Parameters.Experiment.ParticleConcentration = ...
    particle_volume_fraction * [1, 1];

% Constants for diffusion 
channel_width_microns = 5E3;
channel_depth_microns = 10;
flow_rate_ul_min = 0.5;
% T_kelvin = 10E3;
T_kelvin = 0;
viscosity_pas = 1.12E-3;
dx_target_pix = sqrt(tx_pix^2 + ty_pix^2);

% Particle diffusion
particle_diffusion_std_dev_pix = calculate_particle_diffusion(...
    channel_width_microns, channel_depth_microns, ...
    flow_rate_ul_min, objective_magnification, ...
    pixel_size_microns, T_kelvin, dp_microns, viscosity_pas, dx_target_pix);

% Temperature in kelvin
DefaultJob.Parameters.Experiment.ChannelWidth = channel_width_microns;
DefaultJob.Parameters.Experiment.ChannelDepth = channel_depth_microns;
DefaultJob.Parameters.Experiment.FlowRate = flow_rate_ul_min;
DefaultJob.Parameters.Experiment.Temperature = T_kelvin;
DefaultJob.Parameters.Experiment.Viscosity = viscosity_pas;

% Case 1
SegmentItem = DefaultJob;
SegmentItem.SetType = 'mc';
SegmentItem.Parameters.Image.Height = 1024;
SegmentItem.Parameters.Image.Width  = 1024;
SegmentItem.Parameters.Noise.Mean = 0 * [0.1, 0.1];
SegmentItem.Parameters.Noise.Std = 0.03 * [1, 1];

% Update the flow rate
 flow_rate_ul_min = 2.0;
 SegmentItem.Parameters.Experiment.FlowRate = flow_rate_ul_min;
 SegmentItem.CaseName = sprintf('piv_test_running_ensmeble_q_%0.1f_ul_min', flow_rate_ul_min);
%  SegmentItem.CaseName = 'shear_test_01';

 % Particle diffusion
 [particle_diffusion_std_dev_pix, fps] = calculate_particle_diffusion(...
     channel_width_microns, channel_depth_microns, ...
     flow_rate_ul_min, objective_magnification, ...
     pixel_size_microns, T_kelvin, dp_microns, viscosity_pas, dx_target_pix);
 
 
%  % Update diffusion
 SegmentItem.Parameters.Experiment.DiffusionStdDev = ...
     particle_diffusion_std_dev_pix * [1, 1];
%  
  % Update diffusion
%  SegmentItem.Parameters.Experiment.DiffusionStdDev = ...
%      0 * [1, 1];
 
 % Append the job
 JOBLIST(1) = SegmentItem;

%%%%%%%%%%

end








