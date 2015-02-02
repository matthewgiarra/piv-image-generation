function compile_generateParticleImage_3D();
% compile_generateParticleImage_3D(); 
% This function compiles the Matlab function generateParticleImage_3D to
% mex code.
%
% Inputs : None
%
% Outputs: None
%
% Usage: compile_generateParticleImage_3D;
%
% See Also:
%   generateParticleImage_3D

% Example variables
image_height = 1024;
image_width = 1024;
image_depth = 1024;
particle_diameters = coder.typeof(1.00, [inf, 1]);
particle_max_intensities = coder.typeof(1.00, [inf, 1]);
particle_positions_rows = coder.typeof(1.00, [inf, 1]);
particle_positions_columns = coder.typeof(1.00, [inf, 1]);
particle_position_depths = coder.typeof(1.00, [inf, 1]);

% Set up the coder configuration
cfg = coder.config('mex');
cfg.DynamicMemoryAllocation = 'AllVariableSizeArrays';
cfg.GenerateReport = true;

% Run coder to generate the mex file.
codegen -config cfg generateParticleImage_3D -args {image_height, image_width, image_depth, particle_positions_columns, particle_positions_rows, particle_position_depths, particle_diameters, particle_max_intensities};

end