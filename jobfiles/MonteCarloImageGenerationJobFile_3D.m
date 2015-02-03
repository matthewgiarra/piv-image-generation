function JOBLIST = MonteCarloImageGenerationJobFile_3D()

% Region dimensions
regionHeight = 64;
regionWidth = 64;
regionDepth = 64;

DefaultJob.JobOptions.NumberOfProcessors = 1;
DefaultJob.JobOptions.NumberOfDigits = 6;
DefaultJob.JobOptions.RotationRangeType = 'lin';
DefaultJob.JobOptions.RotationAngleUnits = 'rad';
DefaultJob.JobOptions.RunCompiled = 1;
DefaultJob.JobOptions.ReSeed = 1;

DefaultJob.ImageType = 'synthetic';
DefaultJob.SetType = 'mc';
DefaultJob.CaseName = '2015-02-02_translation_only';
DefaultJob.ProjectRepository = '~/Desktop/piv_test_images';

DefaultJob.Parameters.RegionHeight = regionHeight;
DefaultJob.Parameters.RegionWidth = regionWidth;
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1;
DefaultJob.Parameters.Sets.ImagesPerSet = 1;

% Rigid-body displacements (pixels)
DefaultJob.Parameters.TX = 1 * regionWidth  / 8 * [-1, 1];
DefaultJob.Parameters.TY = 1 * regionHeight / 8 * [-1, 1];
DefaultJob.Parameters.TZ = 1 * regionDepth  / 8 * [-1, 1];

% Range of isotropic scaling factors
DefaultJob.Parameters.Scaling = [1 1]; 

% Range of rotation angles (degrees)
DefaultJob.Parameters.Rotation = 0 * [1 1];

% Range of particle concentrations (particles per pixel)
DefaultJob.Parameters.ParticleConcentration = 2E-4 * [1 1];

% Particle diameter (pixels)
DefaultJob.Parameters.ParticleDiameter.Mean = sqrt(8);
DefaultJob.Parameters.ParticleDiameter.Std = [0, 0];

% Noise parameters
DefaultJob.Parameters.Noise.Mean = 0.05;
DefaultJob.Parameters.Noise.Std = 0.05;

% Case 1
SegmentItem = DefaultJob;
SegmentItem.SetType = 'lin';
SegmentItem.CaseName = '2015-02-02_volume_test_lin';
SegmentItem.Parameters.RegionHeight = 64;
SegmentItem.Parameters.RegionWidth = 64;
SegmentItem.Parameters.RegionDepth = 64;
SegmentItem.Parameters.TX =  0 * [0 1];
SegmentItem.Parameters.TY =  0 * [0 1];
SegmentItem.Parameters.TZ =  0 * [1 1];
SegmentItem.Parameters.Rotation_Z_01 = pi/30 * [0 0];
SegmentItem.Parameters.Rotation_Y    = pi/30 * [0 0];
SegmentItem.Parameters.Rotation_Z_02 = pi/6 * [0 0];
SegmentItem.Parameters.Scaling = 0.99 * [1, 1];
SegmentItem.Parameters.Noise.Mean = 0.00;
SegmentItem.Parameters.Noise.Std = 0.00;
SegmentItem.Parameters.ParticleDiameter.Mean = sqrt(8);
SegmentItem.parameters.ParticleDiameter.Std = [0.00, 0.00];
JOBLIST(1) = SegmentItem;


end








