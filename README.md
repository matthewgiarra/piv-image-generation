# PIV Image Generation
These are matlab codes that generate artificial images for evaluating particle image velocimetry (PIV) algorithms. 

This documentation requires some cleanup, but I no longer have MATLAB to run the code myself. If anyone has a solution or wants to buy me a license, I'd be happy to do it. I have not tried running it in Octave.

![png_600x600](https://github.com/user-attachments/assets/ba5734a1-34ee-49b4-aa6d-a26bf5f721e9)

# Overview
`piv-image-generation` is a library of Matlab codes used to generate synthetic photographs for testing planar particle image velocimetry codes (PIV). Particle positions are generated pseudo-randomly, and then advected according to known displacement fields (for full-field images) or coordinate transformations (for Monte–Carlo analysis). The code also saves additional Mat files that contain the ground-truth solutions describing the particle motions, vorticity, etc. These solution files are the basis for assessing the performance of PIV codes.  Particle images are rendered using the algorithm of Olsen and Adrian [1], and the specific details of the implementation can be found in Brady et al. [2]. 

The general workflow for using these codes is as follows. 

1. The user-defined parameters for generating images are specified in `jobfiles`, which are Matlab functions that return a structure or list of structures whose fields specify the parameters that define an image generation job. These parameters include the paths at which to save images, the number of images to be generated, the sizes of the images to be generated, etc.
2. The structures returned by `jobfiles` are passed as arguments to functions that generate and save the synthetic images.
3. Running an image generation job produces the specified set of images, as well as ground-truth solution files for assessing the performance of PIV estimation algorithms.

# Installation
Just clone this repository:

```git clone https://github.com/matthewgiarra/piv-image-generation```

# Usage
## Compiling C functions

We include a C routine for rendering particle images (`c_codes/generateParticleImage.c`), which compiles to a Matlab Mex functions, and runs about 100x faster than the corresponding Matlab version (`core/generateParticleImage.m`).

To compile the `C` routine to a Matlab Mex file, run `compile_scripts/compile_generateParticleImage.m` in Matlab.

```Matlab
cd compile_scripts
compile_generateParticleImage();

% Move the mex file to `piv-image-generation/core`
movefile generateParticleImage.mexa64 ../core
```

If `generateParticleImage.mexa64` exists on the Matlab path, then calling `generateParticleImage()` will run it. If the mex file doesn't exists on the Matlab path, `generateParticleImage.m` will be run instead. No changes to the invoking code are necessary.

# Generating images
## Lamb-Oseen Vortex Full-Field Images
This section provides a guided example for generating synthetic PIV images of a [Lamb-Oseen vortex](https://en.wikipedia.org/wiki/Lamb%E2%80%93Oseen_vortex).

In Matlab, from the `piv-image-generation` root directory:

1. Add the jobfiles directory to the Matlab path.
```matlab
addpath jobfiles;
```

2. Open the PIV image generation job file.
```matlab
open lambOseenVortexImageGenerationJobFile;
```

3. Run the jobfile and return its output to a structure variable called JobFile.

```matlab
JobFile = lambOseenVortexImageGenerationJobFile;
```

4. Pass the `JobFile` structure as an input to `generateImages_LambVortex()`

```matlab
generateImages_LambVortex(JobFile);
```

5.	Inspect the images to make sure they were output properly.

# Reading ground-truth solutions
Lamb Vortex Full-Field Images
The ground-truth Eulerian displacement can be calculated at the center of each PIV interrogation region (IR) by running the Lamb vortex velocity function and inputting the coordinates of the IR centers as the initial computational particle positions. The following code serves as a template example for extracting the ground-truth Eulerian displacement and vorticity fields at the centers of PIV interrogation regions.

## Full Example for reading ground truth solution

```Matlab
% Create a List of the image numbers used in each correlation pair.
firstImageNumbers = startImage : frameStep : endImage;
secondImageNumbers = firstImageNumbers + correlationStep;

% Specify the path to the image generation parameters file
% that was saved concurrently with the synthetic images 
parametersPath = '/path/to/image/generation/parameters/file.mat';

% Load the image generation parameters file.
imageGenerationParameters = load(parametersPath);

% Extract the vortex parameters 
% from the image generation parameters file.
vortexParameters = imageGenerationParameters.VortexParameters;

% Extract the image times from the parameters file.
% Image times in sort-of physical units (i.e. "seconds")
imageTimes = imageGenerationParameters.T;

% Determine the solution times corresponding to the frame numbers
firstImageTimes  = imageTimes(firstImageNumbers);
secondImageTimes = imageTimes(secondImageNumbers);

% Simulated time elapsed between the two images
interFrameTime = secondImageTimes(1) - firstImageTimes(1);

% Set the solution time to halfway between the image times.
% This results in a second-order accurate estimate
% of the Eulerian displacement.
solutionTimes = firstImageTimes + (secondImageTimes - firstImageTimes) / 2;

% This is the number of solution times calculated.
number_of_solution_times = length(solutionTimes);

% These are column vectors corresponding to the
% row (Y) and column (X) coordinates
% of the centers of the PIV interrogation regions
% The variables X and Y are expected to come from
% the external PIV software's solution file.
gridPointsX = X(:);
gridPointsY = Y(:);

% regionWidth and regionHeight are the width and height
% of each interrogation region in pixels. In this example,
% the IRs will be assumed to be 64x64.
regionWidth  = 64;
regionHeight = 64;

% These are the coordinates of the geometric centroids of each PIV window.
% This corresponds to the physical location in the flow that the PIV
% correlation is supposed to measure. For odd sized windows, the geometric
% centroid of the window is at the center pixel. For even-sized windows, it
% is 0.5 pixels to the right of the pixel located at (regionHeight/2) or
% (regionWidth/2). 
xCenter_01 = gridPointsX + 0.5 * (1 - mod(regionWidth,  2));
yCenter_01 = gridPointsY + 0.5 * (1 - mod(regionHeight, 2));

% Vector containing all the x and y grid points.
% This is an input to the vortex velocity function,
% and this format (a single column) is required by ODE45.
gridPointsVector = cat(1, gridPointsX, gridPointsY);

% These lines allocate matrices for the analytical
% velocity fields calculated at each solution time.
%
% Horizontal velocity component
uTrue = zeros([size(X), number_of_solution_times]);

% Vertical velocity component
vTrue = zeros([size(X), number_of_solution_times]);

% Out of plane component of vorticity
rTrue = zeros([size(X), number_of_solution_times]);

% This loops over the different solution times
% and calculates the analytical displacement at each time.
for t = 1 : number_of_solution_times
	
	% This calculates the true Eulerian velocity field for the t'th field.
	[trueVelocities, trueVorticity] = ...
	lambOseenVortexRingVelocityFunction(solutionTimes(t), ...
	[xCenter_01; yCenter_01], vortexParameters);
			
	% This is the ground-truth horizontal component
	% of the velocity field expressed as a vector.
	uTrueVect = interFrameTime * trueVelocities(1 : length(trueVelocities) / 2);
	
	% This is the ground-truth vertical component
	% of the velocity field expressed as a vector.
	vTrueVect = interFrameTime * trueVelocities(length(trueVelocities)/2 + 1 : end);
	
	% This is the ground-truth out-of-plane component
	% of the vorticity field expressed as a vector.
	rTrueVect = interFrameTime * trueVorticity(:);
	
	% This reshapes the horizontal component
	% of the ground-truth velocity field into
	% a matrix of the same size as the 
	% input coordinates.
	uTrue(:, :, t) = flipud(reshape(uTrueVect, size(X)));
	
	% This reshapes the vertical component
	% of the ground-truth velocity field into
	% a matrix of the same size as the 
	% input coordinates.
	vTrue(:, :, t) = -1 * flipud(reshape(vTrueVect, size(X)));
	
	% This reshapes the out-of-plane component
	% of the ground-truth vorticity field into
	% a matrix of the same size as the 
	% input coordinates.
	rTrue(:, :, t) = flipud(reshape(rTrueVect, size(X)));
	
end
```
