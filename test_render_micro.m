

% Parse structure: image dimensions in pixels
Parameters.Image.Width = 1280;
Parameters.Image.Height = 1024;

% Image magnification (microns per pixel)
Parameters.Image.PixelSize = 20;

% Optics parameters
%
% Objective lens
% Magnification (unitless)
Parameters.Optics.Objective.Magnification = 50;

% Objective lens focal length in microns
Parameters.Optics.Objective.FocalLength = 4E3;

% Objective lens numerical aperture (unitless)
Parameters.Optics.Objective.NA = 0.75;

% Laser wavelength in microns
Parameters.Optics.Laser.Wavelength = 0.532;

% Fraction above background intensity to render particles
Parameters.Experiment.IntensityFraction = 0.0;

% Image transformation
Parameters.Transformation = eye(3,3);

% Experiment parameters
% 
% Channel depth in microns
Parameters.Experiment.ChannelDepth = 50;

% Particle diameter in microns
Parameters.Experiment.ParticleDiameter = 1.0;

% Particle concentration (particles per µm^3)
Parameters.Experiment.ParticleConcentration = 5E-3;

% Diffusion
Parameters.Experiment.DiffusionStDev = 0;

generateImagePair_micro_mc(Parameters, 0);