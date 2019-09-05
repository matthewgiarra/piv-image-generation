function compile_all();
% This function compiles all the compile-able codes in this library
% to MEX files. It just calls each of the individual functions' 
% compilation scripts sequentially.
% INPUTS
%   none
%
% OUTPUTS
%   none

% Inform the user of progress
fprintf(['\nCompiling function compile_generateParticleImage.m '...
    'to mex file...\n\n']);
% Compile the 2D Monte Carlo image pair generation code.
compile_generateParticleImage;

% Inform the user of progress
fprintf(['\n\nCompiling function compile_generateParticleImage_3D.m '...
    'to mex file...\n\n']);
% Compile the 3D Monte Carlo image pair generation code.
compile_generateParticleImage_3D;

end