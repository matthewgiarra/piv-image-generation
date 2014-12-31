function JOBLIST = MonteCarloImageGenerationJobFile()

% Region dimensions
regionHeight = 64;
regionWidth = 64;

DefaultJob.JobOptions.NumberOfProcessors = 1;
DefaultJob.JobOptions.NumberOfDigits = 6;
DefaultJob.JobOptions.RotationRangeType = 'lin';
DefaultJob.JobOptions.RotationAngleUnits = 'rad';
DefaultJob.JobOptions.RunCompiled = 1;

DefaultJob.ImageType = 'synthetic';
DefaultJob.SetType = 'mc';
DefaultJob.CaseName = '2014-11-21_translation_only';
DefaultJob.ProjectRepository = '~/Desktop/piv_test_images';

DefaultJob.Parameters.RegionHeight = regionHeight;
DefaultJob.Parameters.RegionWidth = regionWidth;
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1;
DefaultJob.Parameters.Sets.ImagesPerSet = 100;

% Rigid-body displacements (pixels)
DefaultJob.Parameters.TX =  1 * regionWidth / 8  * [-1 1];
DefaultJob.Parameters.TY =  1 * regionHeight / 8 * [-1 1];

% Range of isotropic scaling factors
DefaultJob.Parameters.Scaling = [1 1]; 

% Range of rotation angles (degrees)
DefaultJob.Parameters.Rotation = 0 * [1 1];

% Range of Shear rates (pixels per pixel)
DefaultJob.Parameters.ShearX = [0 0]; % Range of horizontal shear
DefaultJob.Parameters.ShearY = [0 0];  % Range of vertical shear

% Range of particle concentrations (particles per pixel)
DefaultJob.Parameters.ParticleConcentration = [0.025 0.025];

% Particle diameter (pixels)
DefaultJob.Parameters.ParticleDiameter.Mean = sqrt(8);
DefaultJob.Parameters.ParticleDiameter.Std = [0, 0];

% Noise parameters
DefaultJob.Parameters.Noise.Mean = 0.05;
DefaultJob.Parameters.Noise.Std = 0.05;

% Case 1
SegmentItem = DefaultJob;
SegmentItem.SetType = 'lin';
SegmentItem.CaseName = '2014-12-04_disparity_test_lin';
SegmentItem.Parameters.RegionHeight = 64;
SegmentItem.Parameters.RegionWidth = 64;
SegmentItem.Parameters.TX =  8 * [0 1];
SegmentItem.Parameters.TY =  8 * [0 1];
SegmentItem.Parameters.Rotation = pi/6 * [0 1];
SegmentItem.Parameters.Scaling = [1, 1];
SegmentItem.Parameters.Noise.Mean = 0.05;
SegmentItem.Parameters.Noise.Std = 0.05;
SegmentItem.Parameters.ShearX = 0.00 * [-1 1];
SegmentItem.Parameters.ShearY = 0.00 * [-1 1];
SegmentItem.Parameters.ParticleDiameter.Mean = sqrt(8);
SegmentItem.parameters.ParticleDiameter.Std = [0.00, 0.00];
JOBLIST(1) = SegmentItem;


end








