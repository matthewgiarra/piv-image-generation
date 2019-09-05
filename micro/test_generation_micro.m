
addpath ~/Desktop/spectral-phase-correlation/scripts
addpath ~/Desktop/spectral-phase-correlation/jobfiles
addpath ~/Desktop/spectral-phase-correlation/filtering
addpath ~/Desktop/piv-image-generation
addpath ~/Desktop/piv-image-generation/jobfiles
addpath ~/Desktop/FreezeColors

% Microscope objective
objective_magnification = 60;

% Region sizes
region_height = 128;
region_width = 128;

% Number of images
num_images = 1;

% Displacements
sx = 15.00 + rand;
sy = 0 * (15.00 + 1 * 0.2345);

% Diffusion
diffusion_stdev = 2.1;

% Mean particle diameter                               
d_mean = 0.1;

% Particle concentration
% (volume of particles / volume of domain)
c = 5E-5;

% Image noise
% image_noise_std = 0.03;
image_noise_std = 0.03;

% Gaussian window fraction
g_fract = 0.5;

% Channel depth in microns
channel_depth_microns = 100;

% Window
g = gaussianWindowFilter([region_height, region_width], g_fract * [1, 1], 'fraction');

% No window
% g = ones(region_height, region_width);

% Make the joblist
image_gen_joblist = MonteCarloImageGenerationJobFile_micro;

% Update the image generation joblist
image_gen_joblist(1).Parameters.Image.Height = region_height;
image_gen_joblist(1).Parameters.Image.Width = region_width;
image_gen_joblist(1).Parameters.Experiment.ParticleDiameter = d_mean * [1, 1];
image_gen_joblist(1).Parameters.Experiment.DiffusionStdDev = diffusion_stdev * [1, 1];
image_gen_joblist(1).Parameters.Translation.X = sx * [1, 1];
image_gen_joblist(1).Parameters.Translation.Y = sy * [1, 1];
image_gen_joblist(1).Parameters.Sets.ImagesPerSet = num_images;
image_gen_joblist(1).Parameters.Experiment.ParticleConcentration = c * [1, 1];
image_gen_joblist(1).Parameters.Noise.Std = image_noise_std * [1, 1];
image_gen_joblist(1).Parameters.Optics.Objective.Magnification = objective_magnification;
image_gen_joblist(1).Parameters.Experiment.ChannelDepth = channel_depth_microns;

% % SPC job list
spc_joblist = spcJobList_mc;

% Update file names
spc_joblist.CaseName = image_gen_joblist(1).CaseName;
spc_joblist.Parameters.RegionHeight = region_height;
spc_joblist.Parameters.RegionWidth = region_width;

% Update the SPC joblist
spc_joblist.Parameters.Images.End = num_images;
spc_joblist.Parameters.Sets.ImagesPerSet = num_images;
spc_joblist.Parameters.Processing.EnsembleLength = num_images;

% Generate the images
[imageMatrix1, imageMatrix2] = generateMonteCarloImageSet_micro(image_gen_joblist(1));
% load('~/Desktop/piv_test_images/analysis/data/synthetic/mc/piv_test_constant_diffusion/128x128/raw/mc_h128_w128_00001/raw/raw_image_matrix_mc_h128_w128_seg_000001_000010.mat');

% Run the SPC job list

% Window the images
image_01 = double(imageMatrix1(:, :, 1));
image_02 = double(imageMatrix2(:, :, 1));

% FFTs
f1 = fft2((image_01 - mean(image_01(:))) .* g);
f2 = fft2((image_02 - mean(image_02(:))) .* g);

% Correlate
cc = f1 .* conj(f2);

% Phase only filter
pc = phaseOnlyFilter(cc);

% Energy filter
spectral_energy_filter = spectralEnergyFilter(region_height, region_width, sqrt(8));

% RPC-filtered plane
rpc = fftshift(abs(real(ifft2(pc .* fftshift(spectral_energy_filter)))));

scc = fftshift(abs(real(ifft2(cc))));
gcc = fftshift(abs(real(ifft2(pc))));

figure(1)
surf(rpc ./ max(rpc(:)));
title('RPC plane of single images');

% Do the analysis
job_save_path_list = runMonteCarloCorrelationJobFile(spc_joblist);

% Load the results
load(job_save_path_list{1}{1});

% Calculate the errors
tx_err_scc = (TX_TRUE - tx_scc);
ty_err_scc = (TY_TRUE - ty_scc);
err_mag_scc = sqrt(ty_err_scc.^2 + tx_err_scc.^2);

tx_err_rpc = (TX_TRUE - tx_rpc);
ty_err_rpc = (TY_TRUE - ty_rpc);
err_mag_rpc = sqrt(ty_err_rpc.^2 + tx_err_rpc.^2);

tx_err_apc = (TX_TRUE - tx_apc);
ty_err_apc = (TY_TRUE - ty_apc);
err_mag_apc = sqrt(ty_err_apc.^2 + tx_err_apc.^2);

figure(2);
h1 = subplot(1, 2, 1);
imagesc(image_01); axis image
axis off
title(sprintf('$t_x = %0.2f, t_y = %0.2f$', sx, sy), ...
    'interpreter', 'latex', ...
    'FontSize', 16);
caxis([0, intmax('uint8')]);
colormap gray;

h2 = subplot(1, 2, 2);
imagesc(fftshift(angle(pc)));
axis image
axis off
title({'Phase Correlation', sprintf('Diffusion = %0.2f', diffusion_stdev)}, ...
    'interpreter', 'latex', ...
    'FontSize', 16);

p = get(gca, 'position');
p(1) = 0.52;
set(gca, 'position', p);


% correlation_plot_save_name = sprintf('ensemble_correlation_plot_h_%d_w_%d_tx_%0.2f_ty_%0.2f_diff_%0.2f.png', ...
%     region_height, region_width, sx, sy, diffusion_stdev);
% 
% correlation_plot_save_dir = '~/Desktop/ensemble_plots/correlation_plots';
% 
% if ~exist(correlation_plot_save_dir, 'dir')
%     mkdir(correlation_plot_save_dir);
% end

% correlation_plot_save_path = fullfile(correlation_plot_save_dir, correlation_plot_save_name);
% print(1, '-dpng', '-r300', correlation_plot_save_path);


figure(3);
plot(err_mag_scc, '-k', 'linewidth', 2);
hold on
plot(err_mag_rpc, '-r', 'linewidth', 2);
plot(err_mag_apc, '-b', 'linewidth', 2);
hold off
h = legend('SCC', 'RPC', 'APC'); 
set(h, 'fontsize', 16);
set(gca, 'FontSize', 16);
axis square
xlabel('Number of pairs', 'FontSize', 16);
ylabel('Translation error magnitude (pixels)', 'FontSize', 16);
title({['$t_x = ' num2str(sx, '%0.2f') ', t_y = ' ...
    num2str(sy, '%0.2f') '\, \textrm{pix}$'], ['$ \left(u^{\prime}, v^{\prime}\right) = ' ...
    num2str(diffusion_stdev, '%0.2f') ' \, \textrm{pix/frame} $'], ...
    ['$ \bar{d_p} = ' ...
    num2str(d_mean, '%0.2f') ' \, \textrm{pix} $']
    } , 'interpreter', 'latex');

ylim([0, 1]);

% error_plot_save_name = sprintf('ensemble_plot_h_%d_w_%d_tx_%0.2f_ty_%0.2f_diff_%0.2f.png', ...
%     region_height, region_width, sx, sy, diffusion_stdev); 
% 
% error_plot_save_dir = '~/Desktop/ensemble_plots/error_plots';
% 
% if ~exist(error_plot_save_dir, 'dir')
%     mkdir(error_plot_save_dir);
% end
% % 
% % error_plot_save_path = fullfile(error_plot_save_dir, error_plot_save_name);
% % 
% % print(2, '-dpng', '-r300', error_plot_save_path);
% 
% 












