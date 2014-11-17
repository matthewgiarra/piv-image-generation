function generateMonteCarloImageSet(JOBLIST)
% generateMonteCarloImageSet(SIZE, STARTSET, ENDSET, NIMAGES, NPROCESSORS)
% This code generates sets of pairs of particle pattern images that are
% deformed by affine transformations.
%
% INPUTS
%   SIZE = 2 x 1 element vector of integers specifying the height and width (in pixels) of the
%   images to be generated, i.e., SIZE = [HEIGHT WIDTH]. HEIGHT and WIDTH
%   much each be positive integers greater than zero. 
%
%   STARTSET = Integer specifying the number of the  first set of images to
%   generate. 
%
%   ENDSET = Integer specifying the number of the final set of images to
%   generate.
%
%   NIMAGES = Integer specifying the number of image pairs to generate in
%   each set of images.
%
%   NPROCESSORS = Integer specifying the number of processor cores to use for
%    image generation. 
%
% OUTPUTS
%   This code does not save any outputs to the workspace. Instead, it saves
%   images to locations specified within the body of the code.
%
% SEE ALSO
%   generateImagePair.m

% This is the beginning of the user inputs.

% Number of jobs in the job list
nJobs = length(JOBLIST);

% Loop over all the jobs.
for n = 1 : nJobs
    
    % Read the jobfile from the job list
    JobFile = JOBLIST(n);

    % Image Repository. This is line specifies where the image folders will be saved. 
    projectRepository = JobFile.ProjectRepository;
    imageType = JobFile.ImageType;
    setType = JobFile.SetType;
    caseName = JobFile.CaseName;

    % Start / end sets
    startSet = JobFile.Parameters.Sets.Start;
    endSet = JobFile.Parameters.Sets.End;
    imagesPerSet = JobFile.Parameters.Sets.ImagesPerSet;

    % Number of digits in the image set names.
    numberOfDigits = JobFile.JobOptions.NumberOfDigits;
    numberFormat = ['%0' num2str(numberOfDigits) '.0f'];

    % Number of processors
    nProcessors = JobFile.JobOptions.NumberOfProcessors;
    
    % Flag for whether or not to run compiled codes
    run_compiled = JobFile.JobOptions.RunCompiled;

    % Height and width of subregion images
    regionHeight = JobFile.Parameters.RegionHeight;
    regionWidth = JobFile.Parameters.RegionWidth;

    % Path to the image save directory
    imageSaveDir = fullfile(projectRepository, 'analysis', 'data', ...
        imageType, setType, caseName, [num2str(regionHeight) 'x' num2str(regionWidth)], 'raw');
    
    % Create the repository directory if it doesn't exist already
    if ~exist(imageSaveDir, 'dir')
        mkdir(imageSaveDir)
    end

    % Rigid-body displacements (pixels)
    tX = JobFile.Parameters.TX;
    tY = JobFile.Parameters.TY;

    % Range of isotropic scaling factors
    scaling = JobFile.Parameters.Scaling; 

    % Range of logs-10 of rotation angles (degrees)
    rotation_angle_range = JobFile.Parameters.Rotation;
    rotation_range_type = JobFile.JobOptions.RotationRangeType;
    rotation_angle_units = JobFile.JobOptions.RotationAngleUnits;

    % Range of Shear rates (pixels per pixel)
    shearX = JobFile.Parameters.ShearX; % Range of horizontal shear
    shearY = JobFile.Parameters.ShearY;  % Range of vertical shear

    % Range of particle concentrations (particles per pixel)
    concentration = JobFile.Parameters.ParticleConcentration;

    % Particle diameter (pixels)
    particle_diameter_std  = JobFile.Parameters.ParticleDiameter.Std;
    particle_diameter_mean = JobFile.Parameters.ParticleDiameter.Mean;
    
    % Noise parameters
    noiseMean = JobFile.Parameters.Noise.Mean;
    noiseStd  = JobFile.Parameters.Noise.Std;

    % Specify explicitly the bounds of the scaling parameter for the Monte
    % Carlo simulation
    Si = scaling(1);
    Sf = scaling(2);

    % Specify explicitly the bounds of the rotation angles
    Ri = rotation_angle_range(1);
    Rf = rotation_angle_range(2);

    % Specify the type of rotation range (linear or logarithmic)
    if regexpi(rotation_range_type, 'lin');
        rotationRangeType = 'lin';
        isLin = 1;
        isLog = 0;
    elseif regexpi(rotation_range_type, 'log');
        rotationRangeType = 'log';
        isLin = 0;
        isLog = 1;
    end

    % Specify the units for the rotation angles (degrees or radians)
    if regexpi(rotation_angle_units, 'deg')
        rotationAngleUnits = 'deg';
        isDeg = 1;
        isRad = 0;
    elseif regexpi(rotation_angle_units, 'rad')
        rotationAngleUnits = 'rad';
        isDeg = 0;
        isRad = 1;
    end

    % Specify explicitly the bounds of the horizontal shearing parameter
    ShearXi = shearX(1);
    ShearXf = shearX(2);

    % Specify explicitly the bounds of the vertical shearing parameter
    ShearYi = shearY(1);
    ShearYf = shearY(2);

    % Specify explicitly the bounds of the horizontal displacement parameter
    TXi = tX(1);
    TXf = tX(2);

    % Specify explicitly the bounds of the vertical displacement parameter
    TYi = tY(1);
    TYf = tY(2);

    % Specify explicitly the bounds of the particle concentration parameter
    Conci = concentration(1);
    Concf = concentration(2);
    
    % Specify explicitly the bounds of the particle diameter standard
    % deviation 
    particle_diameter_std_i = particle_diameter_std(1);
    particle_diameter_std_f = particle_diameter_std(2);

    % Specify the Image type
    % mc means monte carlo analysis
    % h and w are the image heights and widths in pixels
    casePrefix = [setType '_h' num2str(regionHeight) '_w' num2str(regionWidth)];

    % Number of image sets
    nSets = endSet - startSet + 1;

    % Monte carlo or linear progress
    isMC = ~isempty(regexpi(setType, 'mc'));
    isLinProg = ~isempty(regexpi(setType, 'lin'));

    % Generate the images. Loop over all of the specified image sets.
    for s = 1:nSets
        disp(['Generating set ' num2str(s) ' of ' num2str(nSets)]);
        fprintf(1, 'Image Set %04.0f\n', s); % Inform the user by printing a message to the screen.
        imageDir = fullfile(imageSaveDir, [casePrefix '_' num2str(startSet + s - 1, '%05.0f')], 'raw'); % Directory in which to save images
        parameterDir = fullfile(imageSaveDir, [casePrefix '_' num2str(startSet + s - 1, '%05.0f')], 'parameters'); % Directory in which to save parameters file

        % Make raw image directory if it doesn't already exist
        if ~exist(imageDir, 'dir')
            mkdir(imageDir);
        end

        % Make parameter directory if it doesn't already exist
        if ~exist(parameterDir, 'dir')
            mkdir(parameterDir);
        end

        % Specify the file path to the Parameters file that will be saved.     
        parametersFilePath = fullfile(parameterDir, ['imageParameters_' setType '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_seg_' num2str(1, numberFormat) '_' num2str(imagesPerSet, numberFormat) '.mat']);

        % Specify the file path to the image matrix.
        imageMatrixFilePath = fullfile(imageDir, ['raw_image_matrix_' setType '_h' num2str(regionHeight) '_w' num2str(regionWidth) '_seg_' num2str(1, numberFormat) '_' num2str(imagesPerSet, numberFormat) '.mat'] );

        % Particle diameters
        % Change this to a variable size
        
        % Linear distribution image sets.
        if isLinProg
            % Make linear distribution of transformation parameters.
            Parameters.Concentration = (linspace(Conci, Concf, imagesPerSet))';
            Parameters.Scaling = (linspace(Si, Sf, imagesPerSet))';
            rotationAnglesRaw = (linspace(Ri, Rf, imagesPerSet))';
            Parameters.ShearX = (linspace(ShearXi, ShearXf, imagesPerSet))';
            Parameters.ShearY = (linspace(ShearYi, ShearYf, imagesPerSet))';
            Parameters.TranslationX = (linspace(TXi, TXf, imagesPerSet))';
            Parameters.TranslationY = (linspace(TYi, TYf, imagesPerSet))';  
            Parameters.ParticleDiameterStd = (linspace(particle_diameter_std_i, particle_diameter_std_f, imagesPerSet))';

        % Monte carlo image sets
        elseif isMC
            % Initialize all of the fields in the parameters array
            Parameters.Concentration = zeros(imagesPerSet, 1);
            Parameters.Scaling = zeros(imagesPerSet, 1);
            rotationAnglesRaw = zeros(imagesPerSet, 1);
            Parameters.Rotation = zeros(imagesPerSet, 1);
            Parameters.ShearX = zeros(imagesPerSet, 1);
            Parameters.ShearY = zeros(imagesPerSet, 1);
            Parameters.TranslationX = zeros(imagesPerSet, 1);
            Parameters.TranslationY = zeros(imagesPerSet, 1);
            Parameters.ParticleDiameterStd = zeros(imagesPerSet, 1);
            Parameters.Tforms = zeros(3, 3, imagesPerSet);

            % Populate the transformation parameters
            for k = 1 : imagesPerSet

                Parameters.Concentration(k) = Conci + (Concf - Conci) * rand; % Random concentration
                Parameters.Scaling(k) =  Si + (Sf - Si) * rand; % Random scaling;
                rotationAnglesRaw(k) = Ri + (Rf - Ri) * rand;

                Parameters.ShearX(k) = ShearXi + (ShearXf - ShearXi) * rand; % Random horizontal shearing;
                Parameters.ShearY(k) = ShearYi + (ShearYf - ShearYi) * rand; % Random vertical shearing;
                Parameters.TranslationX(k) = TXi + (TXf - TXi) * rand; % Random horizontal displacement;
                Parameters.TranslationY(k) =  TYi + (TYf - TYi) * rand; % Random vertical displacement;
                Parameters.ParticleDiameterStd(k) = particle_diameter_std_i + (particle_diameter_std_f - particle_diameter_std_i) * rand;
            end
        end

        % Modify Rotation angles to fit the specified units and type
        if isLin % Case of linear range
            rotationAngles = rotationAnglesRaw;
        elseif isLog % Case of logarithmic range
            rotationAngles = 10 .^ rotationAnglesRaw;
        end

        % Set the rotation angle in the vector of angles
        if isRad % Case of radians
            Parameters.Rotation = rotationAngles;
        elseif isDeg % Case of degrees
            Parameters.Rotation = deg2rad(rotationAngles);
        end

        % Create the image transformations
        for k = 1 : imagesPerSet
            Parameters.Tforms(:, :, k) = makeAffineTransform(Parameters.Scaling(k), Parameters.Rotation(k),...
            Parameters.ShearX(k), Parameters.ShearY(k), Parameters.TranslationX(k), Parameters.TranslationY(k));     % Generate the affine transformation matrix 
        end

        % Save image height and width to structure
        Parameters.ImageHeight = regionHeight;
        Parameters.ImageWidth = regionWidth;
        Parameters.ScalingRange = scaling;
        Parameters.RotationRange = rotation_angle_range;
        Parameters.ShearXRange = shearX;
        Parameters.ShearYRange = shearY;
        Parameters.TranslationXRange = tX;
        Parameters.TranslationYRange = tY;
        Parameters.ConcentrationRange = concentration;
        Parameters.RotationRangeType = rotationRangeType;
        Parameters.RotationAngleUnits = rotationAngleUnits;
        Parameters.ParticleDiameterStdRange = particle_diameter_std;
        Parameters.ParticleDiameterMean = particle_diameter_mean;
        
        % Save parameters to their own variables to cut down on data transfer with the parallel for loop
        concentrations = Parameters.Concentration;
        particle_diameter_std_list = Parameters.ParticleDiameterStd;
        tforms = Parameters.Tforms;
        
        % Don't specify particle diameters here.
        % Instead draw them from a normal distribution
        % whose standard deviation is drawn from a uniform
        % distribution of possible values.
%         particleDiameters = Parameters.ParticleDiameter;
        
        % Preallocate memory for the image matrices.
        imageMatrix1 = zeros(regionHeight, regionWidth, imagesPerSet, 'uint16');
        imageMatrix2 = zeros(regionHeight, regionWidth, imagesPerSet, 'uint16');

        % Max value of the images
        maxVal = double(intmax(class(imageMatrix1)));

        % Make noise matrices. The 2.8 corresponds to the multiple of the standard
        % deviation corresponding to a 99.5% coverage factor.
        noiseMatrix1 = noiseMean * maxVal + noiseStd / 2.8 * maxVal * randn(size(imageMatrix1));
        noiseMatrix2 = noiseMean * maxVal + noiseStd / 2.8 * maxVal * randn(size(imageMatrix2));

        % In the case of parallel processing ...
        if nProcessors > 1

            % Start a timer.
            a = tic; 

            % Parallel processing
            parfor k = 1:imagesPerSet % Parallel loop over all the images

                % Generate the image pairs. 
                % Run compiled image generation code
                if run_compiled
                % Run compiled image generation code
                    [image_01, image_02] = generateImagePair_mc_mex(regionHeight, regionWidth, particle_diameter_mean, particle_diameter_std_list(k), concentrations(k), tforms(:, :, k));
                else
                 % Run image generation code
                    [image_01, image_02] = generateImagePair_mc(regionHeight, regionWidth, particle_diameter_mean, particle_diameter_std_list(k), concentrations(k), tforms(:, :, k));
                end

                % Save images to data matrix.;
                imageMatrix1(:, :, k) = image_01 + cast(noiseMatrix1(:, :, k), 'like', image_01);
                imageMatrix2(:, :, k) = image_02 + cast(noiseMatrix2(:, :, k), 'like', image_02);

            end

            % Stop timer
            toc(a);  

        else    

            % Start a timer
            a = tic; 

            % Single thread procesing
            for k = 1:imagesPerSet % Parallel loop over all the images

                % Generate the image pairs. 
                % Run compiled image generation code
                if run_compiled
                    [image_01, image_02] = generateImagePair_mc_mex(regionHeight, regionWidth, particle_diameter_mean, particle_diameter_std_list(k), concentrations(k), tforms(:, :, k));
                else
                 % Run image generation code
                    [image_01, image_02] = generateImagePair_mc(regionHeight, regionWidth, particle_diameter_mean, particle_diameter_std_list(k), concentrations(k), tforms(:, :, k));
                end

                % Save images to data matrix.;
                imageMatrix1(:, :, k) = image_01 + cast(noiseMatrix1(:, :, k), 'like', image_01);
                imageMatrix2(:, :, k) = image_02 + cast(noiseMatrix2(:, :, k), 'like', image_02);

            end

            % Stop timer
            toc(a);  
        end

        % Save the parameters array.
        save(parametersFilePath, 'Parameters');

        % Save the image matrices.
        save(imageMatrixFilePath, 'imageMatrix1', 'imageMatrix2');

        % Inform the user of the save path
        disp(['Saved images to ' imageMatrixFilePath]); 

    end % End ( for s = 1 : nSets )

end % End (for n = 1 : nJobs )

end

