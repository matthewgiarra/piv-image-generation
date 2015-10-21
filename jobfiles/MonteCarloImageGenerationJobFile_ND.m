function JOBLIST = MonteCarloImageGenerationJobFile_ND()

% Region dimensions
regionHeight = 256;
regionWidth  = 256;
regionDepth  = 64;

DefaultJob.JobOptions.ParallelProcessing = 1;
DefaultJob.JobOptions.NumberOfDigits = 6;
DefaultJob.JobOptions.RotationRangeType = 'lin';
DefaultJob.JobOptions.RotationAngleUnits = 'rad';
DefaultJob.JobOptions.RunCompiled = 1;
DefaultJob.JobOptions.ReSeed = 1;

DefaultJob.ImageType = 'synthetic';
DefaultJob.SetType = 'mc';
DefaultJob.CaseName = '2015-04-07_uncertainty_calibration_images';
DefaultJob.ProjectRepository = '~/Desktop/uncertainty_calibration_images';

DefaultJob.Parameters.RegionHeight = regionHeight;
DefaultJob.Parameters.RegionWidth = regionWidth;
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1;
DefaultJob.Parameters.Sets.ImagesPerSet = 10;

% Rigid-body displacements (pixels)
DefaultJob.Parameters.TX = 1 * regionWidth  / 8 * [-1, 1];
DefaultJob.Parameters.TY = 1 * regionHeight / 8 * [-1, 1];
DefaultJob.Parameters.TZ = 1 * regionDepth  / 8 * [-1, 1];

% Range of isotropic scaling factors
DefaultJob.Parameters.Scaling = [1 1]; 

% Range of rotation angles (degrees)
DefaultJob.Parameters.Rotation_Z_01 = 0 * [1 1];
DefaultJob.Parameters.Rotation_Y = 0 * [1 1];
DefaultJob.Parameters.Rotation_Z_02 = 0 * [1 1];

% Range of horizontal shearing factors
DefaultJob.Parameters.ShearX = 1 * [1, 1];

% Range of vertical shearing factors
DefaultJob.Parameters.ShearY = 1 * [1, 1];

% Range of particle concentrations (particles per pixel)
DefaultJob.Parameters.ParticleConcentration = 1000 * [1 1];

% Particle diameter (pixels)
DefaultJob.Parameters.ParticleDiameter.Mean = sqrt(8);
DefaultJob.Parameters.ParticleDiameter.Std = [0, 0];

% Noise parameters
DefaultJob.Parameters.Noise.Mean = 0.05;
DefaultJob.Parameters.Noise.Std = 0.05;

% Image dimensionality
DefaultJob.Parameters.ImageDimensionalitiy = 2;

% Standard deviation of the simulated Gaussian beam profile
% in the Z-direction (world coordinates).
DefaultJob.Parameters.BeamPlaneStdDev = 1E-2;

% Camera parameters
DefaultJob.Parameters.Cameras = default_camera_parameters();

% Domain
DefaultJob.Parameters.Domain.X = 1 * [-1, 1];
DefaultJob.Parameters.Domain.Y = 1 * [-1, 1];
DefaultJob.Parameters.Domain.Z = 0.0 * [0.9, 1.1];

% Case 1
SegmentItem = DefaultJob;
SegmentItem.SetType = 'mc';
SegmentItem.CaseName = '2015-10-20_pinhole_test';
SegmentItem.Parameters.RegionHeight = 1024;
SegmentItem.Parameters.RegionWidth  = 1280;
SegmentItem.Parameters.RegionDepth = 256;
SegmentItem.Parameters.TX =  0.00 * [1 1];
SegmentItem.Parameters.TY =  0 * [1 1];
SegmentItem.Parameters.TZ =  0.1 * [1 1];
SegmentItem.Parameters.Rotation_Z_01 = 0 * pi * [0 1];
SegmentItem.Parameters.Rotation_Y    = 0 * pi * [0 1];
SegmentItem.Parameters.Rotation_Z_02 = 0 * pi * [0 1];
SegmentItem.Parameters.ShearX = 0 * [0.01, 0.2];
SegmentItem.Parameters.ShearY = 0 * [0.01, 0.2];
SegmentItem.Parameters.Scaling = 1 * [1, 1];
SegmentItem.Parameters.ImageNoise.Mean = [0.00, 0.25];
SegmentItem.Parameters.ImageNoise.StdDev = [0.00, 0.10];
SegmentItem.Parameters.ParticleDiameter.Mean = sqrt(8) * [1, 1];
SegmentItem.Parameters.ParticleDiameter.Std = 0  * [1, 1];
SegmentItem.Parameters.ParticleConcentration = 1E4 * [1, 1];
JOBLIST(1) = SegmentItem;

end








