function compile_generateParticleImage_m();
% compile_generateParticleImage_3D(); 
% This function compiles the Matlab function generateParticleImage_3D to
% mex code.
%
% Inputs : None
%
% Outputs: None
%
% Usage: compile_generateParticleImage;
%
% See Also:
%   generateParticleImage.m

% Example variables
image_height = 1024;
image_width = 1024;
particle_diameters = coder.typeof(1.00, [inf, 1]);
particle_max_intensities = coder.typeof(1.00, [inf, 1]);
particle_positions_rows = coder.typeof(1.00, [inf, 1]);
particle_positions_columns = coder.typeof(1.00, [inf, 1]);

% Set up the coder configuration
cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.GenerateReport = true;

% Run coder to generate the mex file.
codegen -config cfg generateParticleImage -args {image_height, image_width, particle_positions_columns, particle_positions_rows, particle_diameters, particle_max_intensities};

end