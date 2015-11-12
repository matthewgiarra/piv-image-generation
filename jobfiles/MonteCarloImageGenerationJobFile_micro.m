function JOBLIST = MonteCarloImageGenerationJobFile_micro()

% Region dimensions
region_width_pixels  = 1280;
region_height_pixels  = 1024;

DefaultJob.JobOptions.ParallelProcessing = 0;
DefaultJob.JobOptions.NumberOfDigits = 6;
DefaultJob.JobOptions.ReSeed = 1;

DefaultJob.ImageType = 'synthetic';
DefaultJob.SetType = 'mc';
DefaultJob.CaseName = 'piv_test_images';
DefaultJob.ProjectRepository = '~/Desktop/piv_test_images';

DefaultJob.Parameters.Image.Height = region_height_pixels;
DefaultJob.Parameters.Image.Width = region_width_pixels;
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1;
DefaultJob.Parameters.Sets.ImagesPerSet = 8;

% Rigid-body displacements (pixels)
DefaultJob.Parameters.Translation.X =  10 * [1 1];
DefaultJob.Parameters.Translation.Y =  0 * [1 1];
DefaultJob.Parameters.Translation.Z =  0 * [0 2];

% Range of isotropic scaling factors
DefaultJob.Parameters.Scaling = 1 * [1 1]; 

% Range of rotation angles (degrees)
DefaultJob.Parameters.Rotation.Z_01 = 0 * [1 1];
DefaultJob.Parameters.Rotation.Y    = 0 * [1 1];
DefaultJob.Parameters.Rotation.Z_02 = 0 * [1 1];

% Sharing
DefaultJob.Parameters.Shear.XY = 0 * [0.00, 0.10];
DefaultJob.Parameters.Shear.XZ = 0.5 * [1, 1];
DefaultJob.Parameters.Shear.YX = 0 * [0.00, 0.10];
DefaultJob.Parameters.Shear.YZ = 0 * [0.00, 0.10];
DefaultJob.Parameters.Shear.ZX = 0 * [0.00, 0.10];
DefaultJob.Parameters.Shear.ZY = 0 * [0.00, 0.10];

% Range of particle concentrations (particles per pixel)
DefaultJob.Parameters.ParticleConcentration = 5E-3 * [1 1];

% Noise parameters
DefaultJob.Parameters.Noise.Mean = 0.05;
DefaultJob.Parameters.Noise.Std = 0.05;

% Optics
% Image magnification (microns per pixel)
DefaultJob.Parameters.Image.PixelSize = 20;

% Optics parameters
% Objective lens
% Magnification (unitless)
DefaultJob.Parameters.Optics.Objective.Magnification = 50;

% Objective lens focal length in microns
DefaultJob.Parameters.Optics.Objective.FocalLength = 4E3;

% Objective lens numerical aperture (unitless)
DefaultJob.Parameters.Optics.Objective.NA = 0.75;

% Laser wavelength in microns
DefaultJob.Parameters.Optics.Laser.Wavelength = 0.532;

% Fraction above background intensity to render particles
DefaultJob.Parameters.Experiment.IntensityFraction = 0.0;
 
% Channel depth in microns
DefaultJob.Parameters.Experiment.ChannelDepth = 50;

% Particle diameter in microns
DefaultJob.Parameters.Experiment.ParticleDiameter = 1.0 * [1, 1];

% Particle concentration (particles per µm^3)
DefaultJob.Parameters.Experiment.ParticleConcentration = 5E-3;

% Diffusion
DefaultJob.Parameters.Experiment.DiffusionStdDev = 3 * [0, 0];

% Case 1
SegmentItem = DefaultJob;
SegmentItem.SetType = 'mc';
SegmentItem.CaseName = 'piv_test_micro';
SegmentItem.Parameters.Image.Height = 512;
SegmentItem.Parameters.Image.Width  = 512;
SegmentItem.Parameters.Scaling = 1 * [1, 1];
SegmentItem.Parameters.ImageNoise.Mean = 0 * [0.1, 0.1];
SegmentItem.Parameters.ImageNoise.StdDev = 0 * [0.10, 0.10];
SegmentItem.Parameters.Experiment.ParticleDiameter = 1 * [1, 1];
SegmentItem.Parameters.ParticleConcentration = 5E-3 * [1, 1];
JOBLIST(1) = SegmentItem;

end








