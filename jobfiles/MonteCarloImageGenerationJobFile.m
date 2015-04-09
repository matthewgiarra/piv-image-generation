function JOBLIST = MonteCarloImageGenerationJobFile()

% Region dimensions
regionHeight = 64;
regionWidth = 64;

% Number of processors to use
% This number doesn't actually 
DefaultJob.JobOptions.ParallelProcessing = 1;

% Number of digits in the file names to save.
DefaultJob.JobOptions.NumberOfDigits = 6;

% Leave this alone for now.
DefaultJob.JobOptions.RotationRangeType = 'lin';

% Units for specifying rotation angles.
DefaultJob.JobOptions.RotationAngleUnits = 'rad';
DefaultJob.JobOptions.RunCompiled = 1;

% Leave as "synthetic"
DefaultJob.ImageType = 'synthetic';

% This specifies whether to sweep through parameters
% linearly (set to 'lin') or to vary them randomly
% in the Monte Carlo sense, choosing values
% from uniform distributions between the ranges
% specified (set to 'mc');
DefaultJob.SetType = 'mc';
DefaultJob.CaseName = '2015-03-16_spc_test';

% This is the path to the parent folder under which 
% the images and parameters will be saved
DefaultJob.ProjectRepository = '~/Desktop/spc_test';

% Height and width of the regions
DefaultJob.Parameters.RegionHeight = regionHeight;
DefaultJob.Parameters.RegionWidth = regionWidth;

% Set numbers
DefaultJob.Parameters.Sets.Start = 1;
DefaultJob.Parameters.Sets.End = 1;
DefaultJob.Parameters.Sets.ImagesPerSet = 10;

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
DefaultJob.Parameters.Noise.Intensity.Mean = 0.05;
DefaultJob.Parameters.Noise.Intensity.Std = 0.05;

% Particle motion noise parameters
% i.e., diffusion.
% These are the range of standard deviations
% of zero-mean Gaussian random noise added to the displacements
% of particles.
% Units are pixels.
DefaultJob.Parameters.DiffusionStdDev = [0 0];

% Case 1
SegmentItem = DefaultJob;
SegmentItem.SetType = 'mc';
SegmentItem.CaseName = '2015-03-13_spc_test';
SegmentItem.Parameters.RegionHeight = 64;
SegmentItem.Parameters.RegionWidth = 64;
SegmentItem.Parameters.TX =  1 * [-5 5];
SegmentItem.Parameters.TY =  1 * [-5 5];
SegmentItem.Parameters.Rotation = pi/6 * [0 0];
SegmentItem.Parameters.Scaling = [1, 1];

% Image noise specified as a fraciton of the full field intensity
% I know this is a bad way to specify this. Sorry.
SegmentItem.Parameters.Noise.Intensity.Mean = 0.1;
SegmentItem.Parameters.Noise.Intensity.Std = 0.1;

SegmentItem.Parameters.ShearX = 0.00 * [-1 1];
SegmentItem.Parameters.ShearY = 0.00 * [-1 1];
SegmentItem.Parameters.ParticleDiameter.Mean = sqrt(8);
SegmentItem.Parameters.ParticleDiameter.Std = [0.00, 0.00];

% Standard deviation of random particle displacements.
SegmentItem.Parameters.DiffusionStdDev = 0 * [2, 2];

% This is the first job in the list of jobs
JOBLIST(1) = SegmentItem;

end








