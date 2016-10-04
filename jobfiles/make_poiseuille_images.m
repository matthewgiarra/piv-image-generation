function make_poiseuille_images(dx_rand_list, num_pairs, img_repo)
% clear;

% Image directory
% img_repo = '~/Desktop/piv_test_images/synthetic';


% Give this case a name
case_name = 'poiseuille';

% Image extension
image_ext = '.tiff';


% % Image stuff
% Number of images
% num_pairs = 10;

% Image dimensions
image_width = 2048;
image_height = 2048;

% Particle positions buffer
x_buffer = -100;
y_buffer = -100;
% Particle position max and min
x_min = 1 + x_buffer;
x_max = image_width  - x_buffer;
y_min = 1 + y_buffer;
y_max = image_height - y_buffer;

% Image noise
noise_mean_fract = 5E-2;
noise_std_fract  = 5E-2;

% Particle stuff
dp_mean = 3;
dp_std =  1;
particle_concentration = 2E-2;

% Magnification in microns per pixel
mag_um_pix = 7.5;

% Light sheet thickness in microns = 4 * std dev
sheet_thickness_microns = 0.5E3;

% Light sheet thickness in pixels
sheet_thickness_pixels = sheet_thickness_microns / mag_um_pix;

% Light sheet standard deviation in pixels
sheet_std_dev_pix = sheet_thickness_pixels / 4;

% Constrain the z coordinates
z_min = -1 * sheet_thickness_pixels / 2;
z_max =      sheet_thickness_pixels / 2;

% Mean displacements
dx_mean = 10 * 2/3;
dy_mean = 0;
dz_mean = 0;

% % list of random displacement standard deviations.
% dx_rand_list = [0.1, 1, 3, 6];
% dx_rand_list = 0.1;

% Number of diffusion cases.
num_diffusion_cases = length(dx_rand_list);

% Remove any dots in the image extension.
image_ext = lower(strrep(image_ext, '.', ''));

% This is a structure of all the data to save
JobFile = {};
JobFile.Parameters.Displacement.Mean.X = dx_mean;
JobFile.Parameters.Displacement.Mean.Y = dy_mean;
JobFile.Parameters.Displacement.Mean.Z = dz_mean;

JobFile.Parameters.Light.SheetThicknessPixels = sheet_thickness_pixels;
JobFile.Parameters.Noise.MeanFraction = noise_mean_fract;
JobFile.Parameters.Noise.StdFraction = noise_std_fract;

JobFile.Parameters.Particles.Concentration = particle_concentration;
JobFile.Parameters.Particles.Diameter.Mean = dp_mean;
JobFile.Parameters.Particles.Diameter.Std = dp_std;

JobFile.Parameters.Image.Height = image_height;
JobFile.Parameters.Image.Width = image_width;


for n = 1 : num_diffusion_cases

    
    
% % Do some werk
% Make lists of images
image_nums_01 = 1 : 2 : 2 * num_pairs;
image_nums_02 = image_nums_01 + 1;

% Random displacement standard deviations
dx_rand_std = dx_rand_list(n);
dy_rand_std = dx_rand_std;
dz_rand_std = dx_rand_std;

% Save random displacements to the job file.
JobFile.Parameters.Displacement.Rand.X = dx_rand_std;
JobFile.Parameters.Displacement.Rand.Y = dy_rand_std;
JobFile.Parameters.Displacement.Rand.Z = dz_rand_std;

% Case name for this diffusion
case_name_current = sprintf('%s_diffusion_%0.2f', case_name, dx_rand_std);

JobFile.CaseName = case_name_current;

% Image directory
image_dir = fullfile(img_repo, case_name_current, 'raw');

% Jobfile directory
job_file_dir = fullfile(img_repo, case_name_current, 'jobfiles');

% Make image directory 
if ~exist(image_dir, 'dir')
    mkdir(image_dir)
end;

% Make jobfile directory 
if ~exist(job_file_dir, 'dir')
    mkdir(job_file_dir)
end;

% Loop over images
parfor k = 1 : num_pairs
    
    % Inform the user
    fprintf(1, 'On image %d of %d\n', k, num_pairs);
    
    % Image save name
    image_name_01 = sprintf('%s_%06d.%s', case_name_current, image_nums_01(k), image_ext);
    image_name_02 = sprintf('%s_%06d.%s', case_name_current, image_nums_02(k), image_ext);
    
    % Image save path
    image_path_01 = fullfile(image_dir, image_name_01);
    image_path_02 = fullfile(image_dir, image_name_02);
    
    % % Make the images
    % Compute the total number of particles
    aug_height = y_max - y_min + 1;
    aug_width = x_max - y_min + 1;
    num_particles = round(particle_concentration * aug_height * aug_width);

    % Particle positions (image 1)
    x_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;
    z_01 = (z_max - z_min) * rand(num_particles, 1) + z_min;

    % Velocities
    dy_rand = dy_rand_std * randn(num_particles, 1);
    dx_rand = dx_rand_std * randn(num_particles, 1);
    dz_rand = dz_rand_std * randn(num_particles, 1);
    
    % Center of the image in the height direction
    yc = image_height / 2;
    
    % Radial coordinate in the height direction
    r = abs((y_01 - yc) / (image_height/2));

    % Total displacements
    dx = - 3/2 * dx_mean * (r.^2 - 1) + dx_rand;
    dy = dy_mean + dy_rand;
    dz = dz_mean + dz_rand;
    
    % Particle positions (image 2)
    x_02 = x_01 + dx;
    y_02 = y_01 + dy;
    z_02 = z_01 + dz;

    % Particle diameters
    dp       = dp_mean + dp_std * randn(num_particles, 1);

    % Particle intensities for the correlated images
    particle_intensities_01 = exp(-z_01.^2 / (2 * sheet_std_dev_pix^2));
    particle_intensities_02 = exp(-z_02.^2 / (2 * sheet_std_dev_pix^2));
    
    % Make the images
    image_01_raw = generateParticleImage(image_height, image_width,...
        x_01, y_01, dp, particle_intensities_01);

    % Generate the second image
    image_02_raw = generateParticleImage(image_height, image_width,...
          x_02, y_02, dp, particle_intensities_02);

    % Noise absolute values
    noise_std = noise_std_fract   * double(intmax('uint8'));
    noise_mean = noise_mean_fract * double(intmax('uint8'));

    % Noise matrices
    noise_mat_01 = noise_mean + noise_std * ...
        randn(image_height, image_width);
    noise_mat_02 = noise_mean + noise_std * ...
        randn(image_height, image_width);

    % Add noise to the images.
    image_01 = image_01_raw;
    image_02 = image_02_raw;
    
    % Calculate the max value of the images
    img_max = max([image_01(:); image_02(:)]);
    
    % Min value of the images
    img_min = min([image_01(:); image_02(:)]);
    
    image_01_shift = image_01 + img_min;
    image_01_scaled = image_01_shift ./ img_max * 0.95 * double(intmax('uint8')) + noise_mat_01;
    
    image_02_shift = image_02 + img_min; % Intentionally done 
    image_02_scaled = image_02_shift ./ img_max * 0.95 * double(intmax('uint8')) + noise_mat_02;
    
    % Save the images
    imwrite(uint8(image_01_scaled), image_path_01, 'compression', 'none');
    imwrite(uint8(image_02_scaled), image_path_02, 'compression', 'none');
    
end

job_file_name = sprintf('%s.mat', case_name_current);
job_file_path = fullfile(job_file_dir, job_file_name);
save(job_file_path, 'JobFile');

end

end

    

