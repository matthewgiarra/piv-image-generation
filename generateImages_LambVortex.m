function  generateImages_LambVortex(JOBFILE)

% generateImagePair_RankineVortex(IMAGEHEIGHT, IMAGEWIDTH, PARTICLEDIAMETER, PARTICLECONCENTRATION, VORTICITY, VORTEXRADIUS, NIMAGES, NPROCESSORS, SAVEDIR, SAVEBASE)
% GENERATEIMAGES(imHeight, imWidth, TRANSFORMATIONMATRIX, saveFlag, saveDirectory)
% Generates a series of synthetic particle images that have been transformed by a
% specified series of transformation matrices. 
%
% The equations used to generate the particles were taken from Equation 16 (on page 5) of the paper
% "Methods for Digital Particle Image Sizing (DPIS): Comparisons and
% Improvements" (2009) by Michael R. Brady, Samuel G. Raben, Pavlos P.
% Vlachos, published in the journal Flow Measurement and Instrumentation.
% (doi:10.1016/j.flowmeasinst.2009.08.001)
%
% INPUTS
%   IMAGEHEIGHT = Height of generated images (pixels)
%
%   IMAGEWIDTH = Width of generate images (pixels)
%
%   PARTICLECONCENTRATION = Average number of particles per pixel
%
%   TRANSFORMATIONMATRIX = 3 x 3 matrix that specifies the affine
%   transformation relating the particle positions in the first image to
%   the particle positions in the second image. 
%
%   SAVEFLAG = Flag to specify whether or not to save the generated images (1 to save, 0 to not save). 
%
%   SAVEDIR = Directory in which to save the generated images. 
%
%   SAVEBASE = Name of the images to save, which will be followed by '_A'
%   and '_B' for the first and second images in the pair. 
%
%   CONTROLCASE = Binary flag (0 or 1) specifying whether to generate a
%   simple six-particle pattern to easily observe the effects of the
%   transformation on the images. This should only be enabled for debugging
%   purposes. 
%
% OUTPUTS
%   This code does not save any outputs to the workspace. Instead, it saves
%   images to locations specified by the inputs to the code.
% 
% SEE ALSO
%    makeAffineTransform

% % % % % % % % % 
% BEGIN FUNCTION %
% % % % % % % % % 

% Load the job file
if ischar(JOBFILE)
    load(JOBFILE)
end

%% Read the job options

% Number of processors
nProcessors = JOBFILE.JobOptions.NumberOfProcessors;

% Flag specifying whether to run compiled codes
run_compiled = JOBFILE.JobOptions.RunCompiled;

% Number of digits in image names
numberOfDigits = JOBFILE.JobOptions.NumberOfDigits;
numberFormat = ['%0' num2str(numberOfDigits) '.0f']; % Numberformat for naming images

% Flag specifying absolute (1) or relative (0) path to image save directory.
projectRepositoryPathIsAbsolute = JOBFILE.JobOptions.DataRepositoryPathIsAbsolute;

% Project repository. This indicates the path to the parent directory
% for the project for which images are being generated.
if projectRepositoryPathIsAbsolute
    projectRepository = JOBFILE.ProjectRepository;
else
    projectRepository = fullfile(determineLocalAetherPath, JOBFILE.ProjectRepository);
end

% Image type. 'experimental' or 'synthetic'.
imageType = JOBFILE.ImageType;

% Image class (uint8 or uint16)
imageClass = JOBFILE.Parameters.ImageClass;

% Set type, i.e., 'vortex'
setType = JOBFILE.SetType;

% Case name, i.e., 'rankinevortex_2013-12-13'
caseName = JOBFILE.CaseName;

% Number of the first set of images to generate
startSet = JOBFILE.Parameters.Sets.Start;

% Number of the last set of images to generate
endSet = JOBFILE.Parameters.Sets.End;

% List of sets
setList = startSet : endSet;

% Number of sets
numberOfSets = length(setList);

% Number of images per set
imagesPerSet = JOBFILE.Parameters.Sets.ImagesPerSet;
 
% Particle concentration (particles per pixel)
particleConcentration = JOBFILE.Parameters.Images.ParticleConcentration;

% Particle diameter (pixels)
particleDiameter = JOBFILE.Parameters.Images.ParticleDiameter; 

% Average particle diameter
particle_diameter_mean = JOBFILE.Parameters.Images.ParticleDiameter.Mean;

% Standard deviation of particle diameter
particle_diameter_standard_deviation = JOBFILE.Parameters.Images.ParticleDiameter.StandardDeviation;

% Height and width of images
imageHeight = JOBFILE.Parameters.Images.Height;
imageWidth = JOBFILE.Parameters.Images.Width;

% Noise mean and standard deviation fractions
noiseMean = JOBFILE.Parameters.Noise.Mean;
noiseStd = JOBFILE.Parameters.Noise.Std;

% Background noise level (fraction of intmax(imageClass))
noiseBackground = JOBFILE.Parameters.Noise.Background;

% Flags for simulating noise and beam
simulateNoise = JOBFILE.JobOptions.SimulateNoise;
simulateBeam = JOBFILE.JobOptions.SimulateBeam;

% Image base name
imageBaseName = ['lambvortex_h' num2str(imageHeight) '_w' num2str(imageWidth) '_'];

% Directory in which to save the images
imageParentDirectory = fullfile(projectRepository, 'analysis', 'data', ...
    imageType, setType, caseName, ['c_' num2str(particleConcentration, '%0.4f')]) ;

% Create image save directory if it doesn't already exist
if ~exist(imageParentDirectory, 'dir')
    mkdir(imageParentDirectory)
end

%%% GENERATE IMAGES

for n = 1 : numberOfSets
    
    % Determine which set is being generated
    currentSet = setList(n);
    
    % Parameters directory
    caseDirectory = fullfile(imageParentDirectory, ['lambvortex_h' num2str(imageHeight) '_w' num2str(imageWidth) '_' num2str(currentSet, '%05.0f')]);
    parametersDir = fullfile(caseDirectory, 'parameters') ;
    imageSaveDir = fullfile(caseDirectory, 'raw');
    
    % Create parameters directory if it doesn't already exist
    if ~exist(parametersDir, 'dir')
        mkdir(parametersDir)
    end
    
    % Create image save directory directory if it doesn't already exist
    if ~exist(imageSaveDir, 'dir')
        mkdir(imageSaveDir)
    end
    
    % Display image save directory
    disp(['Saving images to ' imageSaveDir]);
    
    % Extract vortex parameters
    VortexParameters.CoreRadius = JOBFILE.Parameters.Vortex.CoreRadius;
    VortexParameters.VortexRadius = JOBFILE.Parameters.Vortex.VortexRadius;
    VortexParameters.Angle = JOBFILE.Parameters.Vortex.Angle;
    VortexParameters.PeakVelocity = JOBFILE.Parameters.Vortex.PeakVelocity;
    VortexParameters.PropagationVelocity = JOBFILE.Parameters.Vortex.PropagationVelocity;

    % Time vector
    vortexStartTime = JOBFILE.Parameters.Vortex.StartTime;
    vortexEndTime = JOBFILE.Parameters.Vortex.EndTime;
    vortexTimeVector = linspace(vortexStartTime, vortexEndTime, imagesPerSet);

    % % Image noise level and confidence interval;
    % imageNoiseLevel = JOBFILE.Parameters.Images.ImageNoiseLevel;
    % imageNoiseConfidenceInterval = JOBFILE.Parameters.Images.ImageNoiseConfidenceInterval;

    % Vortex ring center
    VortexParameters.XC = JOBFILE.Parameters.Vortex.Position.X * imageWidth;
    VortexParameters.YC = JOBFILE.Parameters.Vortex.Position.Y * imageHeight;

    % Vortex type.
    VortexParameters.VortexType = 'lamb';

    % Gaussian parameter to control the width of intensity distribution 
    gaussian_light_sheet_standard_deviation = 0.25; 

    % Randomly generate the coordinates of particles
    % Generate horizontal locations of particles in first image (column vector)

    % Left and right extents of the domain over which particles are generated.
    % Particles are generated outside of the image domain so that they can advect outside of the field of view
    % during the particle-path integrations.
    leftParticleLimit = -0.5 * imageWidth;
    rightParticleLimit = 1.5 * imageWidth;

    % Top and bottom extents of the domain over which particles are generated.
    topParticleLimit = -0.5 * imageHeight;
    bottomParticleLimit = 1.5 * imageHeight;

    % Calulate the number of particles to generate
    number_of_particles = round(particleConcentration * (bottomParticleLimit - topParticleLimit) * (rightParticleLimit - leftParticleLimit));

    % Randomly generate the initial positions of particles over the particle domain.
    X0 = leftParticleLimit + (rightParticleLimit - leftParticleLimit) * rand(number_of_particles, 1);
    Y0 = topParticleLimit + (bottomParticleLimit - topParticleLimit) * rand(number_of_particles, 1);

    % Generate lamb vortex particle positions
    disp('Integrating particle positions...');
    [X, Y, T] = lambOseenVortexRing(X0, Y0, vortexTimeVector, VortexParameters);

    % Save the image parameters
    save(fullfile(parametersDir, 'parameters.mat'), 'JOBFILE', 'X', 'Y', 'T', 'imageHeight', 'imageWidth', 'particleDiameter', 'particleConcentration', 'VortexParameters');

    % Uniformly distributed random numbers  from -0.5 to 0.5 corresponding to centers of the Gaussian function
    particleCenters = rand( number_of_particles, 1 ) - 0.5; 

    % Gaussian Function that expresses the intensity distribution of the particle image on the "sensor"
    particleMaxIntensities = (exp( -particleCenters .^ 2 / (2 * gaussian_light_sheet_standard_deviation ^ 2 ) )); 

    % Inform the user.
    disp('Generating particle images...')

    % Beam profile standard deviation, so that FWHM width
    % is equal to half the image width.
    if simulateBeam
        % This is the standard deviation of the in-plane Gaussian 
        % profile of the beam.
        beamStd = imageWidth / (2 * sqrt(2 * log(2)));
        
        % This is the column location of the maxiumum intensity 
        % of the beam.
        beamCenter = imageWidth / 2;
        
        % This creates an array of column locations.
        imageGridX = 1 : imageWidth;
        
        % This calculates the Gaussian profile of the beam intensity
        % as a function of the column position in the image.
        beamProfile = exp(-(imageGridX - beamCenter).^2 / (2 * beamStd.^2));
        
        % This replicates the beam profile in the row direction. 
        % Currently the beam propogation direction is hard-coded to be
        % in the row direction (i.e., top to bottom). 
        beamImage = repmat(beamProfile, [imageHeight, 1]);
    else
        % This sets the beam intensity to be uniform within the image.
        beamImage = ones(imageHeight, imageWidth);
    end

    % Max value of the images
    maxVal = double(intmax(imageClass));

    % noise background image.
    noiseBackgroundImage = noiseBackground .* maxVal .* ones([imageHeight, imageWidth]) .* beamImage;

    % Generate particle diameters
    % Generate the random particle diameters
    particle_diameters = abs(particle_diameter_standard_deviation * randn(number_of_particles, 1) + particle_diameter_mean);
    
    % Generate the images for the subsequent time steps (parallel processing)
    if nProcessors > 1

        % Generate the images.
        parfor k = 1 : imagesPerSet

            % Inform the user
            fprintf(1, ['Generating particle image ' num2str(k) ' of ' num2str(imagesPerSet) '\n']);

            % Generate the subsequent images
            if run_compiled
                ImageOut = generateParticleImage_mex(imageHeight, imageWidth, (X(k, :))', (Y(k, :))', particle_diameters, particleMaxIntensities);
            else
                ImageOut = generateParticleImage(imageHeight, imageWidth, X(k, :), Y(k, :), particle_diameters, particleMaxIntensities);
            end

            % Add noise if requested. Then save the image.
            if simulateNoise
                imwrite(cast((ImageOut ./ max(ImageOut(:)) .* double(intmax(imageClass)) .* beamImage) + (noiseMean * maxVal + noiseStd / 2.8 * maxVal * randn([imageHeight, imageWidth])) + noiseBackgroundImage, imageClass), fullfile(imageSaveDir, [imageBaseName num2str(k, numberFormat) '.tiff']), 'compression', 'none');
            else
            % Save the image
                imwrite(cast(ImageOut ./ max(ImageOut(:)) .* double(intmax(imageClass)) , imageClass ), fullfile(imageSaveDir, [imageBaseName num2str(k, numberFormat) '.tiff']), 'compression', 'none');
            end
        end

    else
    
        % Generate the images for the subsequent time steps (single processor)
        for k = 1 : imagesPerSet 
            % Inform the user
            fprintf(1, ['Generating particle image ' num2str(k) ' of ' num2str(imagesPerSet) '\n']);

            % Generate the subsequent images
            if run_compiled
                ImageOut = generateParticleImage_mex(imageHeight, imageWidth, (X(k, :))', (Y(k, :))', particle_diameters, particleMaxIntensities);
            else
                ImageOut = generateParticleImage(imageHeight, imageWidth, X(k, :), Y(k, :), particle_diameters, particleMaxIntensities);
            end

            % Add noise if requested. Then save the image.
            if simulateNoise
                imwrite(cast((ImageOut ./ max(ImageOut(:)) .* double(intmax(imageClass)) .* beamImage) + (noiseMean * maxVal + noiseStd / 2.8 * maxVal * randn([imageHeight, imageWidth])) + noiseBackgroundImage, imageClass), fullfile(imageSaveDir, [imageBaseName num2str(k, numberFormat) '.tiff']), 'compression', 'none');
            else
            % Save the image
                imwrite(cast(ImageOut ./ max(ImageOut(:)) .* double(intmax(imageClass)) , imageClass ), fullfile(imageSaveDir, [imageBaseName num2str(k, numberFormat) '.tiff']), 'compression', 'none');
            end

        end
    end

end


end % End of function

