function compile_generateImagePair_3D_mc();

% Example variables
IMAGE_HEIGHT = 128;
IMAGE_WIDTH = 128 ;
IMAGE_DEPTH = 128;
PARTICLE_DIAMETER_MEAN = 2.8;
PARTICLE_DIAMETER_STD = 0.00;
PARTICLECONCENTRATION = 0.020;
TRANSFORMATIONMATRIX = eye(4);

% Set up the coder configuration
cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.GenerateReport = true;

% Run coder to generate the mex file.
codegen -config cfg generateImagePair_3D_mc -args {IMAGE_HEIGHT, IMAGE_WIDTH, IMAGE_DEPTH, PARTICLE_DIAMETER_MEAN, PARTICLE_DIAMETER_STD, PARTICLECONCENTRATION, TRANSFORMATIONMATRIX};


end