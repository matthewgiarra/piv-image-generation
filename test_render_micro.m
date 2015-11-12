

% Parse structure: image dimensions in pixels
Parameters.Image.Width = 512;
Parameters.Image.Height = 512;

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

% Transform parameters
R = [0, 0, 0];
S = [1, 1, 1];

% Shearing: [sxy, sxz, syx, syz, szx, szy]
sxy = 0;
sxz = 0;
syx = 0;
syz = 0;
szx = 0;
szy = 0;
SH = [sxy, sxz, syx, syz, szx, szy];
T = [0, 0, 0];
tform = makeAffineTransform_3D(S, R, SH, T);

% Image transformation
Parameters.Transformation = tform;

% Experiment parameters
% 
% Channel depth in microns
Parameters.Experiment.ChannelDepth = 50;

% Particle diameter in microns
Parameters.Experiment.ParticleDiameter = 1.0;

% Particle concentration (particles per µm^3)
Parameters.Experiment.ParticleConcentration = 5E-3;

% Diffusion
Parameters.Experiment.DiffusionStDev = 3;

% Generate image pair.
[IMAGE1, IMAGE2] = generateImagePair_micro_mc(Parameters);

% % Show the images;
% subplot(1, 2, 1);
% imagesc(IMAGE1); axis image; colormap gray;
% c = caxis;
% caxis([0, max(c)]);
% 
% subplot(1, 2, 2);
% imagesc(IMAGE2); axis image; colormap gray;
% caxis([0, max(c)]);


