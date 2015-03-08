function compile_generateImagePair_mc();

% Example variables
IMAGEHEIGHT = 128;
IMAGEWIDTH = 128 ;
PARTICLE_DIAMETER_MEAN = 2.8;
PARTICLE_DIAMETER_STD = 0.00;
PARTICLECONCENTRATION = 0.020;
DIFFUSION_STDEV = 2.00;
TRANSFORMATIONMATRIX = eye(3);

% Set up the coder configuration
cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.GenerateReport = true;

% Run coder to generate the mex file.
codegen -config cfg generateImagePair_mc -args {IMAGEHEIGHT, IMAGEWIDTH, PARTICLE_DIAMETER_MEAN, PARTICLE_DIAMETER_STD, PARTICLECONCENTRATION, DIFFUSION_STDEV, TRANSFORMATIONMATRIX};


end