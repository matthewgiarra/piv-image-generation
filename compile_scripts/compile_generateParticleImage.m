function compile_generateParticleImage();

% Compile the codes using the mac command / linux commands

try
	if isunix && ~ismac % Case for linux machines
	    mex -O CFLAGS="\$CFLAGS -std=c99" generateParticleImage.c
	else
		% This command should work with both mac and windows.
		mex -O CFLAGS="\$CFLAGS -O3" generateParticleImage.c;
	end
	
	% Inform the user
	fprintf(['Compiled generateParticleImage.c to generateParticleImage.' mexext '\n']);

catch
	fprintf('Error compiling codes.\n');
end

end