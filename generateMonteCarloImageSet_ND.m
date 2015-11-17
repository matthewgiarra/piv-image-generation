function [imageMatrix1, imageMatrix2] = generateMonteCarloImageSet_ND(JOBLIST)
% generateMonteCarloImageSet(JOBLIST)
% This code generates sets of pairs of particle pattern images that are
% deformed by affine transformations.
%
% INPUTS
% JOBLIST = Output of the function MonteCarloImageGenerationJobFile_ND.m
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
	
	% World domain
	world_domain_x = JobFile.Parameters.Domain.X;
	world_domain_y = JobFile.Parameters.Domain.Y;
	world_domain_z = JobFile.Parameters.Domain.Z;
	
	% World domain array. The vectorization and transpose operations here
	% ensure that the array is formatted properly, with each 
	% row cooresponding to a different coordinate
	world_domain = [(world_domain_x(:))'; (world_domain_y(:))'; (world_domain_z(:))'];

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
    parallel_processing = JobFile.JobOptions.ParallelProcessing;
	
	% Simulate cameras
	simulate_cameras = JobFile.JobOptions.SimulateCameras;
	
	% Camera parameters
	if simulate_cameras
		
        % Default camera parameters
        JobFile.Parameters.Cameras = default_camera_parameters();
        
        % Update amera parameters
        camera_parameters = calculate_camera_parameters(JobFile.Parameters.Cameras);
        
        % Save to output structure.
        Parameters.Cameras = camera_parameters;
	end
    
    % Flag for whether or not to run compiled codes
    run_compiled = JobFile.JobOptions.RunCompiled;
    
    % Dimensionality of the images (either 2 or 3).
    image_dimensionality = JobFile.Parameters.ImageDimensionalitiy;
    
    % Beam plane standard deviation
    beam_plane_std_dev = JobFile.Parameters.Beam.StdDev;
	beam_plane_center_world_z = JobFile.Parameters.Beam.WorldCenterZ;
    
    % Flag specifying whether or not to re-generate the 
    % coordinates of the particles in the first frame
    % of each pair. 
    %
    % Setting reSeed = 1 generates particle image pairs
    % that are correlated with eachother but not with 
    % other image pairs in the series.
    %
    % Setting reSeed = 0 generates a single series of 
    % particle images that are all correlated with one another,
    % and namely, where the coordinates of particles in the first
    % image are used to calculate the coordinats of particles in all
    % subsequent images. 
    %
    % This option only applies when the option JobFile.SetType == 'lin'.
    % In other words, initial particle positions are automatically
    % re-generated for each frame pair when generating Monte Carlo series,
    % even if reSeed == 1.
    reSeed = JobFile.JobOptions.ReSeed;

    % Height and width of subregion images
    region_height = JobFile.Parameters.RegionHeight;
    region_width  = JobFile.Parameters.RegionWidth;
    region_depth  = JobFile.Parameters.RegionDepth;
    
    % Path to the image save directory
    imageSaveDir = fullfile(projectRepository, 'analysis', 'data', ...
        imageType, setType, caseName, [num2str(region_height) 'x' num2str(region_width)], 'raw');
    
    % Create the repository directory if it doesn't exist already
    if ~exist(imageSaveDir, 'dir')
        mkdir(imageSaveDir)
    end

    % Rigid-body displacements (pixels)
    tX = JobFile.Parameters.TX;
    tY = JobFile.Parameters.TY;
    tZ = JobFile.Parameters.TZ;

    % Range of isotropic scaling factors
    scaling = JobFile.Parameters.Scaling; 
    
    % Range of horizontal shearing
    shearX = JobFile.Parameters.ShearX;
    
    % Range of vertical shearing
    shearY = JobFile.Parameters.ShearY;

    % Ranges of Euler-decomposed rotation angles (degrees)
    rotation_angle_range_Z_01 = JobFile.Parameters.Rotation_Z_01;
    rotation_angle_range_Y    = JobFile.Parameters.Rotation_Y;
    rotation_angle_range_Z_02 = JobFile.Parameters.Rotation_Z_02;
    
    % Type of rotation range (linear or log)
    rotation_range_type = JobFile.JobOptions.RotationRangeType;
    
    % Units of rotation angles (degrees or radians)
    rotation_angle_units = JobFile.JobOptions.RotationAngleUnits;

	% Particle concentration
	concentration = JobFile.Parameters.ParticleConcentration;
    
    % Particle diameter (pixels)
    particle_diameter_std  = JobFile.Parameters.ParticleDiameter.Std;
    particle_diameter_mean = JobFile.Parameters.ParticleDiameter.Mean;
    
    % Noise parameters
    image_noise_mean = JobFile.Parameters.ImageNoise.Mean;
    image_noise_std_dev  = JobFile.Parameters.ImageNoise.StdDev;
	


    % Specify explicitly the bounds of the scaling parameter for the Monte
    % Carlo simulation
    Si = scaling(1);
    Sf = scaling(2);

    % Specify explicitly the bounds of the Euler-decomposed rotation angles
    % First rotation about the Z axis
    Ri_Z_01 = rotation_angle_range_Z_01(1);
    Rf_Z_01 = rotation_angle_range_Z_01(2);
    
    % Rotation about the Y axis
    Ri_Y = rotation_angle_range_Y(1);
    Rf_Y = rotation_angle_range_Y(2);
    
    % Second rotation about the Z axis.
    Ri_Z_02 = rotation_angle_range_Z_02(1);
    Rf_Z_02 = rotation_angle_range_Z_02(2);

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

    % Specify explicitly the bounds of the horizontal shearing parameters
    shear_x_i = shearX(1);
    shear_x_f = shearX(2);
    
    % Specify explicitly the bounds of the vertical shearing parameters
    shear_y_i = shearY(1);
    shear_y_f = shearY(2);
    
    % Specify explicitly the bounds of the horizontal displacement parameter
    TXi = tX(1);
    TXf = tX(2);
   
    % Specify explicitly the bounds of the vertical displacement parameter
    TYi = tY(1);
    TYf = tY(2);
    
    % Specify explicitly the bounds of the depth-wise displacement parameter
    TZi = tZ(1);
    TZf = tZ(2);

    % Specify explicitly the bounds of the particle concentration parameter
    Conci = concentration(1);
    Concf = concentration(2);
	
    % Specify explicitly the bounds of the particle diameter standard
    % deviation 
    particle_diameter_std_i = particle_diameter_std(1);
    particle_diameter_std_f = particle_diameter_std(2);
    
    % Specify the bounds of the particle mean diameter
    particle_diameter_mean_i = particle_diameter_mean(1);
    particle_diameter_mean_f = particle_diameter_mean(2);
    
    % Specify the bounds of the image noise mean intensity
    image_noise_mean_i = image_noise_mean(1);
    image_noise_mean_f = image_noise_mean(2);
    
    % Specify the bounds of the image noise standard deviation
    image_noise_std_dev_i = image_noise_std_dev(1);
    image_noise_std_dev_f = image_noise_std_dev(2);

    % Specify the Image type
    % mc means monte carlo analysis
    % h and w are the image heights and widths in pixels
    casePrefix = [setType '_h' num2str(region_height) ...
                          '_w' num2str(region_width)];

    % Number of image sets
    nSets = endSet - startSet + 1;

    % Monte carlo or linear progress
    isMC = ~isempty(regexpi(setType, 'mc'));
    isLinProg = ~isempty(regexpi(setType, 'lin'));

    % Generate the images. Loop over all of the specified image sets.
    for s = 1 : nSets
        disp(['Generating set ' num2str(s) ' of ' num2str(nSets)]);
        fprintf(1, 'Image Set %04.0f\n', s); % Inform the user by printing a message to the screen.
        
        % Directory in which to save images
        imageDir = fullfile(imageSaveDir, ...
            [casePrefix '_' num2str(startSet + s - 1, '%05.0f')], ...
            'raw'); 
        
        % Directory in which to save parameters file.
        parameterDir = fullfile(imageSaveDir, ...
            [casePrefix '_' num2str(startSet + s - 1, '%05.0f')],...
            'parameters'); 

        % Make raw image directory if it doesn't already exist
        if ~exist(imageDir, 'dir')
            mkdir(imageDir);
        end

        % Make parameter directory if it doesn't already exist
        if ~exist(parameterDir, 'dir')
            mkdir(parameterDir);
        end

        % Specify the file path to the Parameters file that will be saved.     
        parametersFilePath = fullfile(parameterDir, ['imageParameters_' setType '_h' num2str(region_height) '_w' num2str(region_width) '_seg_' num2str(1, numberFormat) '_' num2str(imagesPerSet, numberFormat) '.mat']);

        % Specify the file path to the image matrix.
        imageMatrixFilePath = fullfile(imageDir, ['raw_image_matrix_' setType '_h' num2str(region_height) '_w' num2str(region_width) '_seg_' num2str(1, numberFormat) '_' num2str(imagesPerSet, numberFormat) '.mat'] );
        
        % Linear distribution image sets.
        if isLinProg
            % Make linear distribution of transformation parameters.
            Parameters.Concentration = (linspace(Conci, Concf, imagesPerSet))';
            Parameters.Scaling = (linspace(Si, Sf, imagesPerSet))';
            
            % Euler-decomposed rotation angles
            rotation_angles_Z_01_raw = (linspace(Ri_Z_01, Rf_Z_01, imagesPerSet))';
            rotation_angles_Y_raw =    (linspace(Ri_Y,    Rf_Y,    imagesPerSet))';
            rotation_angles_Z_02_raw = (linspace(Ri_Z_02, Rf_Z_02, imagesPerSet))';
            
            % Concatonate the Euler-decomposed rotation angles
            % into a single array
            rotation_angles_raw = [rotation_angles_Z_01_raw, ...
                               rotation_angles_Y_raw, ...
                               rotation_angles_Z_02_raw];
            
            % Shearing vectors
            Parameters.ShearX = (linspace(shear_x_i, shear_x_f, imagesPerSet))';
            Parameters.ShearY = (linspace(shear_y_i, shear_y_f, imagesPerSet))'; 
            
            % Create particle pattern displacement vectors
            Parameters.TranslationX = (linspace(TXi, TXf, imagesPerSet))';
            Parameters.TranslationY = (linspace(TYi, TYf, imagesPerSet))'; 
            Parameters.TranslationZ = (linspace(TZi, TZf, imagesPerSet))'; 
            Parameters.ParticleDiameterStd = (linspace(particle_diameter_std_i, particle_diameter_std_f, imagesPerSet))';
            Parameters.ParticleDiameterMean = (linspace(particle_diameter_mean_i, particle_diameter_mean_f, imagesPerSet))';
            Parameters.ImageNoise.Mean = (linspace(image_noise_mean_i, image_noise_mean_f, imagesPerSet))';
            Parameters.ImageNoise.StdDev = (linspace(image_noise_std_dev_i, image_noise_std_dev_f, imagesPerSet))';
            
            
        % Monte carlo image sets
        elseif isMC
            % Initialize all of the fields in the parameters array
            Parameters.Concentration = zeros(imagesPerSet, 1);
            Parameters.Scaling = zeros(imagesPerSet, 1);
            Parameters.Rotation = zeros(imagesPerSet, 3);
            Parameters.TranslationX = zeros(imagesPerSet, 1);
            Parameters.TranslationY = zeros(imagesPerSet, 1);
            Parameters.TranslationZ = zeros(imagesPerSet, 1);
            Parameters.ShearX = zeros(imagesPerSet, 1);
            Parameters.ShearY = zeros(imagesPerSet, 1);
            Parameters.ParticleDiameterStd = zeros(imagesPerSet, 1);
            Parameters.ParticleDiameterMean = zeros(imagesPerSet, 1);
            Parameters.Tforms = zeros(4, 4, imagesPerSet);
            Parameters.ImageNoise.Mean = zeros(imagesPerSet, 1);
            Parameters.ImageNoise.StdDev = zeros(imagesPerSet, 1);
            
            % Initialize array for raw rotation angles.
            rotation_angles_raw = zeros(imagesPerSet, 3);

            % Populate the transformation parameters
            for k = 1 : imagesPerSet

                Parameters.Concentration(k) = Conci + (Concf - Conci) * rand; % Random concentration
                Parameters.Scaling(k) =  Si + (Sf - Si) * rand; % Random scaling;
                
                % Euler-decomposed rotation angles
                rotation_angles_Z_01_raw = Ri_Z_01 + (Rf_Z_01 - Ri_Z_01) * rand;
                rotation_angles_Y_raw    = Ri_Y    + (Rf_Y    - Ri_Y)    * rand;
                rotation_angles_Z_02_raw = Ri_Z_02 + (Rf_Z_02 - Ri_Z_02) * rand;
                
                % Concatonate the Euler-decomposed rotation angles
                % into a single array
                rotation_angles_raw(k, :) = [rotation_angles_Z_01_raw,...
                                             rotation_angles_Y_raw,...
                                             rotation_angles_Z_02_raw];

                % Random horizontal displacement;
                Parameters.TranslationX(k) = TXi + (TXf - TXi) * rand;
                
                % Random vertical displacement;
                Parameters.TranslationY(k) =  TYi + (TYf - TYi) * rand; 
                
                % Random depthwise displacement;
                Parameters.TranslationZ(k) =  TZi + (TZf - TZi) * rand; 
                
                % Random horizontal shearing
                Parameters.ShearX(k) = shear_x_i + (shear_x_f - shear_x_i) * rand;
                
                % Random vertical shearing
                Parameters.ShearY(k) = shear_y_i + (shear_y_f - shear_y_i) * rand;
                
                % Random mean particle diameter
                Parameters.ParticleDiameterMean(k) = particle_diameter_mean_i + (particle_diameter_mean_f - particle_diameter_mean_i) * rand;
                
                % Random particle diameter standard deviations.
                Parameters.ParticleDiameterStd(k) = particle_diameter_std_i + (particle_diameter_std_f - particle_diameter_std_i) * rand;
            
                % Intensity noise mean
                Parameters.ImageNoise.Mean(k) = image_noise_mean_i + (image_noise_mean_f - image_noise_mean_i) * rand;
                
                % Intensity noise standard deviation
                Parameters.ImageNoise.StdDev(k) = image_noise_std_dev_i + (image_noise_std_dev_f - image_noise_std_dev_i) * rand;                
            
            end
        end

        % Modify Rotation angles to fit the specified units and type
        if isLin % Case of linear range
            rotation_angles = rotation_angles_raw;
        elseif isLog % Case of logarithmic range
            rotation_angles = 10 .^ rotation_angles_raw;
        end

        % Set the rotation angle in the vector of angles
        if isRad % Case of radians
            Parameters.Rotation = rotation_angles;
        elseif isDeg % Case of degrees
            Parameters.Rotation = deg2rad(rotation_angles);
        end

        % Create the image transformations
        for k = 1 : imagesPerSet
            
            % Generate the 3D similarity transformation matrix
            Parameters.Tforms(:, :, k) = ...
                makeAffineTransform_3D(Parameters.Scaling(k), ...
                                           Parameters.Rotation(k, :), ...
                                           Parameters.ShearX(k), ...
                                           Parameters.ShearY(k), ...
                                           Parameters.TranslationX(k), ...
                                           Parameters.TranslationY(k), ...
                                           Parameters.TranslationZ(k));     
        end

        % Save image height and width to structure
        Parameters.ImageHeight = region_height;
        Parameters.ImageWidth = region_width;
        Parameters.ScalingRange = scaling;
        
        Parameters.TranslationXRange = tX;
        Parameters.TranslationYRange = tY;
        Parameters.TranslationZRange = tZ;
        Parameters.ConcentrationRange = concentration;
        Parameters.RotationRangeType = rotationRangeType;
        Parameters.RotationAngleUnits = rotationAngleUnits;
        Parameters.ParticleDiameterStdRange = particle_diameter_std;
        Parameters.ParticleDiameterMeanRange = particle_diameter_mean;
        Parameters.Beam = JobFile.Parameters.Beam;
		
	   
        % Don't specify particle diameters here.
        % Instead draw them from a normal distribution
        % whose standard deviation is drawn from a uniform
        % distribution of possible values.
        
        % Max value of the images
        % The image class probably shouldn't
        % be hard coded.
        maxVal = double(intmax('uint16'));

        % Save parameters to their own variables to cut down on data transfer with the parallel for loop
        concentrations = Parameters.Concentration;
        particle_diameter_mean_list = Parameters.ParticleDiameterMean;
        
        % Particle diameter standard deviation
        particle_diameter_std_list = Parameters.ParticleDiameterStd;
        tforms = Parameters.Tforms;
        image_noise_mean_list = Parameters.ImageNoise.Mean;
        image_noise_std_dev_list = Parameters.ImageNoise.StdDev .* ...
            image_noise_mean_list;
        
        % Make noise matrix for the series of first images. 
        % The 2.8 corresponds to the multiple of the standard
        % deviation corresponding to a 99.5% coverage factor.
        % This is outside the if-statement because this array
        % gets created whether or not images are re-seeded for each pair.
        %
        % Pick between 2D and 3D
        switch image_dimensionality
            
            % Case of 2D images
            case 2
                
                % Allocate the first noise matrix (2D).
                noiseMatrix1 = zeros(region_height, region_width,...
                    imagesPerSet);
					
                % Allocate the first noise matrix (2D).
                noiseMatrix2 = zeros(region_height, region_width,...
	                    imagesPerSet);
                
                % Populate each array in the noise matrix.
                for k = 1 : imagesPerSet
                
                    % Create the noise matrix for the list of second images.
                    noiseMatrix1(:, :, k) = ...
                        image_noise_mean_list(k) * maxVal + ...
                        image_noise_std_dev_list(k) * maxVal * ...
                        randn([region_height, region_width]);
						
	                % Create the noise matrix for the list of second images.
	                noiseMatrix2(:, :, k) = ...
	                    image_noise_mean_list(k) * maxVal + ...
	                    image_noise_std_dev_list(k) * maxVal * ...
	                    randn([region_height, region_width]);	
						
						
						
                end

            otherwise
                
                % Allocate the noise matrix (3D).
                noiseMatrix1 = zeros(region_height, region_width,...
                    region_depth, imagesPerSet);
					
                % Allocate the noise matrix (3D).
                noiseMatrix2 = zeros(region_height, region_width,...
                    region_depth, imagesPerSet);
                
                % Populate each array in the noise matrix.
                for k = 1 : imagesPerSet

                    % Create the noise matrix for the list of second images.
                    noiseMatrix1(:, :, :, k) = ...
                        image_noise_mean_list(k) * maxVal + ...
                        image_noise_std_dev_list(k) * maxVal * ...
                        randn([region_height, region_width, region_depth]);
						
	                    % Create the noise matrix for the list of second images.
	                    noiseMatrix2(:, :, :, k) = ...
	                        image_noise_mean_list(k) * maxVal + ...
	                        image_noise_std_dev_list(k) * maxVal * ...
	                        randn([region_height, region_width, region_depth]);

                end             
        end
        

        % Start a timer to measure image generation time
        a = tic;
        
        % Generate images, reseeding the initial particle positions
        % for each frame pair.
        if reSeed
        
            % Pick between 2D and 3D
            switch image_dimensionality
                case 2
                    
                    % Preallocate memory for the two image matrices.
                    imageMatrix1_2D = zeros(region_height, region_width, ...
                        imagesPerSet, 'uint16');
                    imageMatrix2_2D = zeros(region_height, region_width, ...
                        imagesPerSet, 'uint16');
                    
                    % Set variables so parfor loop can run
                    imageMatrix1_3D = [];
                    imageMatrix2_3D = [];
 
                otherwise
                    
                    % Preallocate memory for the two image matrices.
                    imageMatrix1_3D = zeros(region_height, region_width, ...
                        region_depth, imagesPerSet, 'uint16');
                    imageMatrix2_3D = zeros(region_height, region_width, ...
                        region_depth, imagesPerSet, 'uint16');
                    
                    % Set variables so parfor loop can run
                    imageMatrix1_2D = [];
                    imageMatrix2_2D = [];
                     
            end
            
            % In the case of parallel processing ...
            if parallel_processing
                % Parallel processing
                parfor k = 1 : imagesPerSet % Parallel loop over all the images
				
				  if simulate_cameras
					% Generate the image pairs.
                      [image_01, image_02] = generateImagePair_2D_pinhole_mc( ...
                          world_domain,  ...
                          particle_diameter_mean_list(k), ...
                          particle_diameter_std_list(k), ...
                          concentrations(k), beam_plane_center_world_z, beam_plane_std_dev, ...
                          tforms(:, :, k), camera_parameters, ...
                          run_compiled);
				  else
	                  % Generate the image pairs.
	                    [image_01, image_02] = generateImagePair_ND_mc( ...
	                        region_height, region_width, region_depth,  ...
	                        particle_diameter_mean_list(k), ...
	                        particle_diameter_std_list(k), ...
	                        concentrations(k), beam_plane_std_dev,...
	                        tforms(:, :, k), image_dimensionality, ...
	                        run_compiled);					  
				  end
						
                    switch image_dimensionality
                        case 2
                                                    
                        % Save z-projection of first volume 
                        % to the data matrix
						imageMatrix1_2D(:, :, k) = ...
                            cast( ...
                                maxVal * mat2gray(sum(image_01, 3)) + ...
                                noiseMatrix1(:, :, k), ...
                                'like', maxVal...
                            );
                        
                         % Save z-projection of second volume
                         % to data matrix
                        imageMatrix2_2D(:, :, k) = ...
                            cast(...
                                maxVal * mat2gray(sum(image_02, 3)) + ...
                                noiseMatrix2(:, :, k), ...
                                'like', maxVal...
                            );
                            
                        otherwise
                            
                            % Save first image to data matrix
                            imageMatrix1_3D(:, :, :, k) = cast(...
                                maxVal * mat2gray(image_01) + ...
                                    noiseMatrix1(:, :, :, k), ...
                                    'like', maxVal...
                                );

                             % Save second image to data matrix
                            imageMatrix2_3D(:, :, :, k) = cast(...
                                maxVal * mat2gray(image_02) + ...
                                    noiseMatrix2(:, :, :, k), ...
                                    'like', maxVal...
                                );
                            
                    end
                    

                end
            else
                  
                % Single thread procesing
                for k = 1 : imagesPerSet % Parallel loop over all the images

  				  if simulate_cameras
  					% Generate the image pairs.
                        [image_01, image_02] = generateImagePair_2D_pinhole_mc( ...
                            world_domain,  ...
                            particle_diameter_mean_list(k), ...
                            particle_diameter_std_list(k), ...
                            concentrations(k), beam_plane_center_world_z, beam_plane_std_dev, ...
                            tforms(:, :, k), camera_parameters, ...
                            run_compiled);
  				  else
  	                  % Generate the image pairs.
  	                    [image_01, image_02] = generateImagePair_ND_mc( ...
  	                        region_height, region_width, region_depth,  ...
  	                        particle_diameter_mean_list(k), ...
  	                        particle_diameter_std_list(k), ...
  	                        concentrations(k), beam_plane_std_dev,...
  	                        tforms(:, :, k), image_dimensionality, ...
  	                        run_compiled);					  
  				  end
                    
                    switch image_dimensionality
                        case 2
                                                    
                        % Save z-projection of first volume 
                        % to the data matrix
                        imageMatrix1_2D(:, :, k) = ...
                            cast( ...
                                maxVal * mat2gray(sum(image_01, 3)) + ...
                                noiseMatrix1(:, :, k), ...
                                'like', maxVal...
                            );
                        
                         % Save z-projection of second volume
                         % to data matrix
                        imageMatrix2_2D(:, :, k) = ...
                            cast(...
                                maxVal * mat2gray(sum(image_02, 3)) + ...
                                noiseMatrix2(:, :, k), ...
                                'like', maxVal...
                            );
                            
                        otherwise
                            
                            % Save first image to data matrix
                            imageMatrix1_3D(:, :, :, k) = cast(...
                                maxVal * mat2gray(image_01) + ...
                                    noiseMatrix1(:, :, :, k), ...
                                    'like', maxVal...
                                );

                             % Save second image to data matrix
                            imageMatrix2_3D(:, :, :, k) = cast(...
                                maxVal * mat2gray(image_02) + ...
                                    noiseMatrix2(:, :, :, k), ...
                                    'like', maxVal...
                                );
                            
                    end
                    
                end 
                

                
            end
            
            % Save the image matrices to the output variables.
            switch image_dimensionality
                case 2
                    imageMatrix1 = imageMatrix1_2D;
                    imageMatrix2 = imageMatrix2_2D;
                otherwise
                    imageMatrix1 = imageMatrix1_3D;
                    imageMatrix2 = imageMatrix2_3D;
            end  

        % Generate images in the linear-progression sense.
        % This statement is accessed when reSeed == 0;
        
        else
            % Update the lists of particle diameter standard deviations
            % and particle concentrations to all have the same value, 
            % since only one was used for each. 
            % Keep these in list format to be compatible with subsequent 
            % codes that are expecting lists rather than a single value. 
            particle_diameter_std_list(:) = particle_diameter_std_list(1);
            concentrations(:) = concentrations(1);
			
			if simulate_cameras
				[imageMatrix1, imageMatrix2] = generateImageSeries_2D_pinhole( ...
                            world_domain,  ...
                            particle_diameter_mean(1), ...
                            particle_diameter_std_list(1), ...
                             concentrations(1), beam_plane_center_world_z, beam_plane_std_dev, ...
                            tforms, camera_parameters, ...
                            run_compiled);
                        
                       
				
			else
	            % Generate the images in the linear-progression sense
	            imageMatrix1 = generateImageSeries_ND(...
	                                      region_height,...
	                                      region_width, ...
	                                      region_depth, ...
	                                      particle_diameter_mean(1),...
	                                      particle_diameter_std_list(1), ...
	                                      concentrations(1), ...
	                                      beam_plane_std_dev, ...
	                                      tforms, image_dimensionality, ...
	                                      parallel_processing, ...
	                                      run_compiled); 
			end

           
           % Add the noise matrix
           imageMatrix1 = uint16(imageMatrix1) + cast(noiseMatrix1(:, :, 1), 'uint16');
           imageMatrix2 = uint16(imageMatrix2) + cast(noiseMatrix2, 'uint16');
           
        end
        
        % Stop timer
        toc(a);
        
        % Save the parameters array.
        save(parametersFilePath, 'Parameters');

        % Save the image matrices.
        save(imageMatrixFilePath, 'imageMatrix1', 'imageMatrix2');

        % Inform the user of the save path
        disp(['Saved images to ' imageMatrixFilePath]); 

    end % End ( for s = 1 : nSets )

end % End (for n = 1 : nJobs )

end

