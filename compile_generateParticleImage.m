function compile_generateParticleImage();

% Example variables
IMAGEHEIGHT = 1024;
IMAGEWIDTH = 1024;
PARTICLEDIAMETER = 2.8;
PARTICLEMAXINTENSITIES = zeros(26000, 1);
X = zeros(26000, 1);
Y = zeros(26000, 1);

% Set up the coder configuration
cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.GenerateReport = true;

% Set the dynamic variable sizes
coder.varsize('X', [inf, 1], [true, false]);
coder.varsize('Y', [inf, 1], [true, false]);
coder.varsize('PARTICLEMAXINTENSITIES', [inf, 1], [true, false]);

% Run coder to generate the mex file.
codegen -config cfg generateParticleImage -args {IMAGEHEIGHT, IMAGEWIDTH, X, Y, PARTICLEDIAMETER, PARTICLEMAXINTENSITIES};

end