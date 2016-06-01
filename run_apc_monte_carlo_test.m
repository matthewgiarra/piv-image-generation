
function run_apc_monte_carlo_test(N, S, generate_images, do_processing)

if nargin < 4
    do_processing = 1;
end

if nargin < 3
    generate_images = 1;
end

if nargin < 2
    S = 1;
end

if nargin < 1
    N = 1;
end

% Get repository paths
spc_repo_src = get_spc_repo('src');
spc_repo_top = get_spc_repo('top');
image_gen_repo_src = get_image_gen_repo('src');

% Add paths SPC
addpath(fullfile(spc_repo_src, 'correlation_algorithms'));
addpath(fullfile(spc_repo_src, 'filtering'));
addpath(fullfile(spc_repo_src, 'jobfiles'));
addpath(fullfile(spc_repo_src, 'phase_unwrapping'));
addpath(fullfile(spc_repo_src, 'scripts'));


% Add image gen paths
addpath(image_gen_repo_src);
addpath(fullfile(image_gen_repo_src, 'jobfiles'));

% Target displacements
tx_pix = 15;
ty_pix = 0;

% Set base name
set_base_name = sprintf('piv_test_running_ensemble_%d_ppf', tx_pix);


% Physical stuff for diffusion
% Constants for diffusion 
channel_width_microns = 5E3;
channel_depth_microns = 100;
T_kelvin = 300;
viscosity_pas = 1.12E-3;
dx_target_pix = sqrt(tx_pix^2 + ty_pix^2);

% Flow rates to test in uL/min
flow_rate_vect = [0.50, 5.0, 50.0, 10.0];

% Load the image generation job list
image_gen_job_list = MonteCarloImageGenerationJobFile_micro();

% Extract the first image gen jobfile
ImageGenJobFile = image_gen_job_list(1);

% Load the SPC job list
piv_job_list = spcJobList_mc();

% Extract the spc job file. Only run the first one.
piv_jobfile = piv_job_list(1);

% Set the PIV jobfile repository path
piv_jobfile.Parameters.RepositoryPath = spc_repo_top;

% Same for the Image gen repository path
ImageGenJobFile.ProjectRepository = spc_repo_top;

% Start and end sets
start_set = 1;
end_set = 100;

% Images per set
images_per_set = 10000;

% Images to analyze
start_image = 1;
end_image = images_per_set;
skip_image = 1;

% Sets vector
set_vect = round(linspace(start_set, end_set + 1, N+1));

% Start and end sets
start_set_current = set_vect(S);
end_set_current = set_vect(S + 1) - 1;

% Sets that the current machine will perform
set_vect_current = start_set_current : end_set_current;

% Sets per jobo on the current machine
num_sets_per_job = length(set_vect_current);

% Number of jobs
num_flow_rates = length(flow_rate_vect);

% Loop over all the jobs
for n = 1 : num_flow_rates
    
    % Current flow rate
    flow_rate_current = flow_rate_vect(n);
   
    % Update the number of images to generate
    ImageGenJobFile.Parameters.Sets.ImagesPerSet = images_per_set;
    
    % Update the number of digits in the file names
    piv_jobfile.JobOptions.NumberOfDigits = ...
        ImageGenJobFile.JobOptions.NumberOfDigits;
    
    % Update the set type
    piv_jobfile.SetType = ImageGenJobFile.SetType;
    piv_jobfile.ImageType = ImageGenJobFile.ImageType;
    
    % Particle diameter in microns
    dp_microns = ImageGenJobFile.Parameters.Experiment.ParticleDiameter(1);
    
    % Objective magnification
    objective_magnification = ...
        ImageGenJobFile.Parameters.Optics.Objective.Magnification;
    
    % Pixel size in microns
    pixel_size_microns = ImageGenJobFile.Parameters.Image.PixelSize;
    
    % Calculate diffusion
    % Particle diffusion
    particle_diffusion_std_dev_pix = calculate_particle_diffusion(...
        channel_width_microns, channel_depth_microns, ...
        flow_rate_vect(n), objective_magnification, ...
        pixel_size_microns, T_kelvin, dp_microns, viscosity_pas, dx_target_pix);
    
    % Set the diffusion in pixels per frame
    ImageGenJobFile.Parameters.Experiment.DiffusionStdDev = ...
        particle_diffusion_std_dev_pix * [1, 1];
    
    % Update the number of sets
    piv_jobfile.Parameters.Sets.ImagesPerSet = images_per_set;
    
    % Update the image numbers
    % Make sure it's starting at 1 and not skipping any.
    piv_jobfile.Parameters.Images.Start = start_image;
    piv_jobfile.Parameters.Images.Skip = skip_image;
    piv_jobfile.Parameters.Images.End = end_image;
    
    % Update the region sizes
    piv_jobfile.Parameters.RegionHeight = ...
        ImageGenJobFile.Parameters.Image.Height;
    piv_jobfile.Parameters.RegionWidth = ...
        ImageGenJobFile.Parameters.Image.Width;
    
    % Case name
    input_case_name = ...
        sprintf('%s_q_%0.1f_ul_min', set_base_name, flow_rate_current);
    
    % Determine the case name of the PIV image gen job
    output_case_name = sprintf('%s_', input_case_name);    
    
    % Set the case names of the jobfiles to be equivalent
    ImageGenJobFile.CaseName = input_case_name;
    piv_jobfile.CaseName = input_case_name;
    
    % Output case name
    piv_jobfile.Filepaths.Output.BaseName = output_case_name;
             			
    % Loop over all the sets
    % Generate a set then run APC
    % rather than generating all the images
    % and then running all the APC
   for s = 1 : num_sets_per_job

        % Copy the set number
        set_num = set_vect_current(s);

        % Copy the image gen jobfile
         I = ImageGenJobFile;

         % Copy the PIV jobfile
         P = piv_jobfile;

        % Generate hte images if requested
        if generate_images

            % Update jobfile with new translations
            I.Parameters.Translation.X = (tx_pix + rand - 0.5) * [1, 1];
            I.Parameters.Translation.Y = (ty_pix + rand - 0.5) * [1, 1];

            % Update set number
            I.Parameters.Sets.Start = set_num;
            I.Parameters.Sets.End = set_num;
            
            % Flow rate
            I.Parameters.Experiment.FlowRate = flow_rate_current;

            % Generate the images
            generateMonteCarloImageSet_micro(I);

        end

        % Process the images if requested
        if do_processing

            % Update the set numbers
            P.Parameters.Sets.Start = set_num;
            P.Parameters.Sets.End = set_num;
            P.OtherInfo.DiffusionStdDev = ...
                I.Parameters.Experiment.DiffusionStdDev(1);
				
            % Run the PIV processing
            runMonteCarloCorrelationJobFile(P);

        end
    end      

end

end


function toplevel = get_top_level();
    if islinux();
        toplevel = '/home/shannon/b/aether';
    else
        toplevel = '/Volumes/aether_b';
    end

end

function repo = get_spc_repo(level)

    if islinux();
        top_repo = ...
    '/home/shannon/b/aether/Projects/APC/';

    switch lower(level)
        case 'top'
            repo = top_repo;
        case 'src'
            repo = fullfile(top_repo, ...
                'analysis', 'src', 'spectral-phase-correlation');
    end

    else
        repo = '~/Desktop/spectral-phase-correlation';
    end

end

function repo = get_image_gen_repo(level)

    if islinux();
        top_repo = ...
    '/home/shannon/b/aether/Projects/piv-image-generation/';

    switch lower(level)
        case 'top'
            repo = top_repo;
        case 'src'
            repo = fullfile(top_repo, 'analysis', 'src', 'piv-image-generation');
    end

    else
        repo = '~/Desktop/piv-image-generation';
    end
	
end












