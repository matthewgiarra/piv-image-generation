function [imageMatrix1, imageMatrix2] = generateMonteCarloImageSet_micro(JOBLIST)
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
	
    % Image Repository. This is line specifies
    % where the image folders will be saved. 
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
    
    % Monte Carlo or linear
    isMC = regexpi(setType, 'mc');
    isLin = regexpi(setType, 'lin');
    
    % Height and width of subregion images
    region_height_pixels   = JobFile.Parameters.Image.Height;
    region_width_pixels    = JobFile.Parameters.Image.Width;
    
    % Depth of the micro channel (into-plane direction)
    channel_depth_microns  = JobFile.Parameters.Experiment.ChannelDepth;
    
    % Pixel size in microns (Assumes square pixels)
    pixel_size_microns = JobFile.Parameters.Image.PixelSize;
    
    % Image bit depth
    image_bit_depth = JobFile.Parameters.Image.BitDepth;
    
    % Objective parameters
    JobFile.Parameters.Optics.Objective = ...
        load_objective_parameters(JobFile.Parameters.Optics.Objective.Name);
    objective_magnification = JobFile.Parameters.Optics.Objective.Magnification;
    focal_length_microns = JobFile.Parameters.Optics.Objective.FocalLength;
    NA = JobFile.Parameters.Optics.Objective.NA;
    working_distance_microns = JobFile.Parameters.Optics.Objective.WorkingDistance;
    
    % Laser stuff
    wavelength_microns = JobFile.Parameters.Optics.Laser.Wavelength;
    
    % Path to the image save directory
    imageSaveDir = fullfile(projectRepository, 'analysis', 'data', ...
        imageType, setType, caseName, ...
        [num2str(region_height_pixels) 'x' ...
        num2str(region_width_pixels)], 'raw');
    
    % Create the repository directory if it doesn't exist already
    if ~exist(imageSaveDir, 'dir')
        mkdir(imageSaveDir)
    end
    
    % Rigid-body displacements (pixels)
    tX = JobFile.Parameters.Translation.X;
    tY = JobFile.Parameters.Translation.Y;
    tZ = JobFile.Parameters.Translation.Z;
    
    % Range of isotropic scaling factors
    scaling = JobFile.Parameters.Scaling;
    
    % Ranges of Euler-decomposed rotation angles (degrees)
    rotation_angle_range_Z_01 = JobFile.Parameters.Rotation.Z_01;
    rotation_angle_range_Y    = JobFile.Parameters.Rotation.Y;
    rotation_angle_range_Z_02 = JobFile.Parameters.Rotation.Z_02;
    
  	% Particle concentration (particles per microliter)
	concentration = JobFile.Parameters.Experiment.ParticleConcentration;
    
    % Particle diameter (pixels)
    particle_diameter_microns  = ...
        JobFile.Parameters.Experiment.ParticleDiameter;
    
    % Particle diffusion
    particle_diffusion_std_dev = JobFile.Parameters.Experiment.DiffusionStdDev;
    
    % Noise parameters
    image_noise_mean     = JobFile.Parameters.Noise.Mean;
    image_noise_std_dev  = JobFile.Parameters.Noise.Std;
    
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
    
    % Specify explicitly the bounds of the horizontal shearing parameters
    shear_xy = JobFile.Parameters.Shear.XY;
    shear_xz = JobFile.Parameters.Shear.XZ;
    shear_yx = JobFile.Parameters.Shear.YX;
    shear_yz = JobFile.Parameters.Shear.YZ;
    shear_zx = JobFile.Parameters.Shear.ZX;
    shear_zy = JobFile.Parameters.Shear.ZY;
    
    % XY shearing limits
    shear_xy_i = shear_xy(1);
    shear_xy_f = shear_xy(2);
    
    % XZ shearing limits
    shear_xz_i = shear_xz(1);
    shear_xz_f = shear_xz(2);
    
    % YX shearing limits
    shear_yx_i = shear_yx(1);
    shear_yx_f = shear_yx(2);
    
    % YZ shearing limits
    shear_yz_i = shear_yz(1);
    shear_yz_f = shear_yz(2);
    
    % ZX shearing limits
    shear_zx_i = shear_zx(1);
    shear_zx_f = shear_zx(2);
    
    % ZY shearing limits
    shear_zy_i = shear_zy(1);
    shear_zy_f = shear_zy(2);

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
    
    % Specify the bounds of the particle mean diameter
    particle_diameter_i = particle_diameter_microns(1);
    particle_diameter_f = particle_diameter_microns(2);
    
    % Specify the bounds of the image noise mean intensity
    image_noise_mean_i = image_noise_mean(1);
    image_noise_mean_f = image_noise_mean(2);
    
    % Specify the bounds of the image noise standard deviation
    image_noise_std_dev_i = image_noise_std_dev(1);
    image_noise_std_dev_f = image_noise_std_dev(2);
    
    % Diffusion
    diffusion_i = particle_diffusion_std_dev(1);
    diffusion_f = particle_diffusion_std_dev(2);
    
    % Specify the Image type
    % mc means monte carlo analysis
    % h and w are the image heights and widths in pixels
    casePrefix = [setType '_h' num2str(region_height_pixels) ...
                          '_w' num2str(region_width_pixels)];
                      
    % Number of image sets
    nSets = endSet - startSet + 1;

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
        parametersFilePath = fullfile(parameterDir, ['imageParameters_' setType '_h' num2str(region_height_pixels) '_w' num2str(region_width_pixels) '_seg_' num2str(1, numberFormat) '_' num2str(imagesPerSet, numberFormat) '.mat']);

        % Specify the file path to the image matrix.
        imageMatrixFilePath = fullfile(imageDir, ['raw_image_matrix_' setType '_h' num2str(region_height_pixels) '_w' num2str(region_width_pixels) '_seg_' num2str(1, numberFormat) '_' num2str(imagesPerSet, numberFormat) '.mat'] );
        
        % Linear distribution image sets.
        if isLin
            % Make linear distribution of transformation parameters.
            Parameters.Concentration = (linspace(Conci, Concf, imagesPerSet))';
            
            % Particle diffusion
            Parameters.Diffusion = (linspace(diffusion_i, ...
                diffusion_f, imagesPerSet))';
           
            % Isotropic scaling
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
                           
            % Save rotations to structure.
            Parameters.Rotation = rotation_angles_raw;
            
            % Shearing vectors
            Parameters.Shear.XY = (linspace(shear_xy_i, shear_xy_f, imagesPerSet))';
            Parameters.Shear.XZ = (linspace(shear_xz_i, shear_xz_f, imagesPerSet))'; 
            Parameters.Shear.YX = (linspace(shear_yx_i, shear_yx_f, imagesPerSet))'; 
            Parameters.Shear.YZ = (linspace(shear_yz_i, shear_yz_f, imagesPerSet))'; 
            Parameters.Shear.ZX = (linspace(shear_zx_i, shear_zx_f, imagesPerSet))'; 
            Parameters.Shear.ZY = (linspace(shear_zy_i, shear_zy_f, imagesPerSet))'; 
            
            % Create particle pattern displacement vectors
            Parameters.Translation.X = (linspace(TXi, TXf, imagesPerSet))';
            Parameters.Translation.Y = (linspace(TYi, TYf, imagesPerSet))'; 
            Parameters.Translation.Z = (linspace(TZi, TZf, imagesPerSet))'; 
            Parameters.ParticleDiameter = (linspace(particle_diameter_i, particle_diameter_f, imagesPerSet))';
            Parameters.ImageNoise.Mean = (linspace(image_noise_mean_i, image_noise_mean_f, imagesPerSet))';
            Parameters.ImageNoise.StdDev = (linspace(image_noise_std_dev_i, image_noise_std_dev_f, imagesPerSet))';
         
        % Monte carlo image sets
        elseif isMC

            % Random particle diffusion
            Parameters.Diffusion = ...
                diffusion_i + (diffusion_f - diffusion_i) * ...
                rand(imagesPerSet, 1);
            
            % Random particle concentration
            Parameters.Concentration = ...
                Conci + (Concf - Conci) * rand(imagesPerSet, 1);
            
            % Isotropic scaling.
            Parameters.Scaling =  ...
                Si + (Sf - Si) * rand(imagesPerSet, 1); % Random scaling;
            
            % Euler-decomposed rotation angles (Z1)
            rotation_angles_Z_01_raw = Ri_Z_01 + ...
                (Rf_Z_01 - Ri_Z_01) * rand(imagesPerSet, 1);
            
             % Euler-decomposed rotation angles (Y)
            rotation_angles_Y_raw    = Ri_Y    + ...
                (Rf_Y    - Ri_Y)    * rand(imagesPerSet, 1);
            
             % Euler-decomposed rotation angles (Z2)
            rotation_angles_Z_02_raw = Ri_Z_02 + ...
                (Rf_Z_02 - Ri_Z_02) * rand(imagesPerSet, 1);
            
            % Combine Euler rotation angles
            rotation_angles_raw = ...
                [rotation_angles_Z_01_raw, ...
                rotation_angles_Y_raw, ...
                rotation_angles_Z_02_raw];
            
            % Save rotations to structure.
            Parameters.Rotation = rotation_angles_raw;
            
            % Random horizontal displacement;
            Parameters.Translation.X = TXi + (TXf - TXi) * ...
                rand(imagesPerSet, 1);
            
            % Random vertical displacement;
            Parameters.Translation.Y =  TYi + (TYf - TYi) * ...
                rand(imagesPerSet, 1);     
            
            % Random depthwise displacement;
            Parameters.Translation.Z =  TZi + (TZf - TZi) * ...
                rand(imagesPerSet, 1); 

            % Random XY shearing
            Parameters.Shear.XY = shear_xy_i + ...
                (shear_xy_f - shear_xy_i) * ...
                rand(imagesPerSet, 1);
            
             % Random XZ shearing
            Parameters.Shear.XZ = shear_xz_i + ...
                (shear_xz_f - shear_xz_i) * ...
                rand(imagesPerSet, 1);
            
            % Random YX shearing
            Parameters.Shear.YX = shear_yx_i + ...
                (shear_yx_f - shear_yx_i) * ...
                rand(imagesPerSet, 1);
            
            % Random YZ shearing
            Parameters.Shear.YZ = shear_yz_i + ...
                (shear_yz_f - shear_yz_i) * ...
                rand(imagesPerSet, 1);
            
            % Random ZX shearing
            Parameters.Shear.ZX = shear_zx_i + ...
                (shear_zx_f - shear_zx_i) * ...
                rand(imagesPerSet, 1);
            
            % Random ZY shearing
            Parameters.Shear.ZY = shear_zy_i + ...
                (shear_zy_f - shear_zy_i) * ...
                rand(imagesPerSet, 1);
         
            % Random mean particle diameter
            Parameters.ParticleDiameter = ...
                particle_diameter_i + ...
                (particle_diameter_f - particle_diameter_i) *...
                rand(imagesPerSet, 1);
            
            % Intensity noise mean
            Parameters.ImageNoise.Mean = ...
                image_noise_mean_i + ...
                (image_noise_mean_f - image_noise_mean_i) *...
                rand(imagesPerSet, 1);

            % Intensity noise standard deviation
            Parameters.ImageNoise.StdDev = ...
                image_noise_std_dev_i + ...
                (image_noise_std_dev_f - image_noise_std_dev_i) * ...
                rand(imagesPerSet, 1);
              
        end

        % Create the image transformations
        for k = 1 : imagesPerSet
                
            % Isotropic scaling
            Scale = Parameters.Scaling(k) * ones(1, 3);
            
          % Rotation
            Rotation = Parameters.Rotation(k, :);
         
            % Shearing
            Shear = [...
                    Parameters.Shear.XY(k), ...
                    Parameters.Shear.XZ(k) ,...
                    Parameters.Shear.YX(k), ...
                    Parameters.Shear.YZ(k), ...
                    Parameters.Shear.ZX(k), ...
                    Parameters.Shear.ZY(k)];
            
            % Translation
            Translation = [...
                          Parameters.Translation.X(k), ...
                          Parameters.Translation.Y(k), ...
                          Parameters.Translation.Z(k)];
                        
            % Make the transform      
             % Generate the 3D similarity transformation matrix
            Parameters.Tforms(:, :, k) = makeAffineTransform_3D(Scale, ...
                                   Rotation, ...
                                   Shear, ...
                                   Translation);
                                          
        end

        % Save image height and width to structure
        Parameters.ImageHeight = region_height_pixels;
        Parameters.ImageWidth = region_width_pixels;
        
        % Save experiment parmaeters 
        Parameters.Experiment = JobFile.Parameters.Experiment;
        Parameters.Optics = JobFile.Parameters.Optics;
        Parameters.Image = JobFile.Parameters.Image;
        
        % Max value of the images
        % The image class probably shouldn't
        % be hard coded.
        maxVal = 2^image_bit_depth - 1;
        
        % Set the data type for saving images
        if image_bit_depth > 8
            
        % Use 16-bit precision for bit depth greater than 8
            data_type_string = 'uint16';
        else
            
        % Use 8-bit precision for bit depth <= 8.
            data_type_string = 'uint8';
        end
      
        % Save parameters to their own variables to cut 
        % down on data transfer with the parallel for loop
        particle_concentrations_list = Parameters.Concentration;
        particle_diameter_microns_list = Parameters.ParticleDiameter;
        
        % Particle diffusion 
        particle_diffusion_list = Parameters.Diffusion;
        
        % Image transformations
        tforms = Parameters.Tforms;
        
        % List of image noise means
        image_noise_mean_list = Parameters.ImageNoise.Mean;
        
        % List of image noise standard deviations
        image_noise_std_dev_list = Parameters.ImageNoise.StdDev;
            
        % Allocate the image matrix (image 1)
        imageMatrix1 = zeros(region_height_pixels, ...
            region_width_pixels, imagesPerSet);
        
        % Allocate the image matrix (image 2)
        imageMatrix2 = zeros(region_height_pixels, ...
            region_width_pixels, imagesPerSet, 'uint16');
         
        % Start a timer
        t1 = tic;
        
        % Populate each array in the noise matrix.
        parfor k = 1 : imagesPerSet

            % Create the noise matrix for the list of second images.
            noiseMatrix1 = ...
            image_noise_mean_list(k) + ...
            image_noise_std_dev_list(k) * ...
            randn([region_height_pixels, region_width_pixels]);

            % Create the noise matrix for the list of second images.
            noiseMatrix2 = ...
            image_noise_mean_list(k) + ...
            image_noise_std_dev_list(k) * ...
            randn([region_height_pixels, region_width_pixels]);	
        
            % Generate the image pair
            % These images are scaled [0, 1]
            [img_01, img_02] = generateImagePair_micro_mc(...
                region_height_pixels, region_width_pixels,...
                particle_diameter_microns_list(k) ,...
                pixel_size_microns, particle_concentrations_list(k), ...
                channel_depth_microns, objective_magnification, NA, ...
                working_distance_microns, focal_length_microns, wavelength_microns, ...
                particle_diffusion_list(k), ...
                tforms(:, :, k));
         
            % Rescale the images and cast as uint16 
            imageMatrix1(:, :, k) = abs(img_01 + noiseMatrix1);
            imageMatrix2(:, :, k) = abs(img_02 + noiseMatrix2);
      
        end
        
        % Scale and re-cast data
        imageMatrix1 = cast(imageMatrix1 .* maxVal, data_type_string);
        imageMatrix2 = cast(imageMatrix2 .* maxVal, data_type_string);
        
        % End the timer
        t2 = toc(t1);
        
        % Calculate seconds per pair
        seconds_per_pair = t2 / imagesPerSet;
        
        % Informt the user
        fprintf('%0.3f seconds for %d pairs\n', t2, imagesPerSet);
        fprintf('%0.3f seconds per pair\n', seconds_per_pair);
        
        % Add the jobfile to the parameters path
        Parameters.JobFile = JobFile;
        
        % Save the parameters array.
        save(parametersFilePath, 'Parameters', '-v7.3');
        
        % Save the image matrices.
        save(imageMatrixFilePath, 'imageMatrix1', 'imageMatrix2', '-v7.3');
        
        % Inform the user of the save path
        disp(['Saved images to ' imageMatrixFilePath]); 
        
    end % End ( for s = 1 : nSets )
    
end % End (for n = 1 : nJobs )

end

