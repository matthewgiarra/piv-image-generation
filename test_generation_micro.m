
addpath ~/Desktop/spectral-phase-correlation/scripts
addpath ~/Desktop/spectral-phase-correlation/jobfiles

% Region sizes
region_height = 128;
region_width = 128;

% Number of images
num_images = 500;

% Displacements
sx = 12 + 1 * rand;
sy = 13 + 1 * rand;

% Diffusion
diffusion_stdev = 10;

% Mean particle diameter
d_mean = 0.1 * sqrt(8);

% Window
g = gaussianWindowFilter([region_height, region_width], [0.5, 0.5], 'fraction');

% No window
% g = ones(region_height, region_width);

% Make the joblist
image_gen_joblist = MonteCarloImageGenerationJobFile_micro;

% SPC job list
spc_joblist = spcJobList_mc;

% Update the image generation joblist
image_gen_joblist(1).Parameters.Image.Height = region_height;
image_gen_joblist(1).Parameters.Image.Width = region_width;
image_gen_joblist(1).Parameters.Experiment.ParticleDiameter = d_mean * [1, 1];
image_gen_joblist(1).Parameters.Experiment.DiffusionStdDev = diffusion_stdev * [1, 1];
image_gen_joblist(1).Parameters.Translation.X = sx * [1, 1];
image_gen_joblist(1).Parameters.Translation.Y = sy * [1, 1];
image_gen_joblist(1).Parameters.Sets.ImagesPerSet = num_images;

% Update the SPC joblist
spc_joblist.Parameters.Images.End = num_images;
spc_joblist.Parameters.Sets.ImagesPerSet = num_images;
spc_joblist.Parameters.Processing.EnsembleLength = num_images;

% Generate the images
[imageMatrix1, imageMatrix2] = generateMonteCarloImageSet_micro(image_gen_joblist(1));
% load('~/Desktop/piv_test_images/analysis/data/synthetic/mc/piv_test_constant_diffusion/128x128/raw/mc_h128_w128_00001/raw/raw_image_matrix_mc_h128_w128_seg_000001_000010.mat');


% Window the images
image_01 = double(imageMatrix1(:, :, 1)) .* g;
image_02 = double(imageMatrix2(:, :, 1)) .* g;

% FFTs
f1 = fft2(image_01);
f2 = fft2(image_02);

% Correlate
cc = f1 .* conj(f2);

% Phase only filter
pc = phaseOnlyFilter(cc);

close all;

figure(1);
subplot(1, 2, 1)
imagesc(image_01); axis image
axis off

subplot(1, 2, 2);
imagesc(fftshift(angle(pc)));
axis image
axis off

figure(2);
% Do the analysis
job_save_path_list = runMonteCarloCorrelationJobFile(spc_joblist);

% Load the results
load(job_save_path_list{1}{1});

% Calculate the errors
tx_err_scc = (TX_TRUE - tx_scc);
ty_err_scc = (-1*TY_TRUE - ty_scc);
err_mag_scc = sqrt(ty_err_scc.^2 + tx_err_scc.^2);

tx_err_rpc = (TX_TRUE - tx_rpc);
ty_err_rpc = (-1*TY_TRUE - ty_rpc);
err_mag_rpc = sqrt(ty_err_rpc.^2 + tx_err_rpc.^2);

tx_err_apc = (TX_TRUE - tx_apc);
ty_err_apc = (-1*TY_TRUE - ty_apc);
err_mag_apc = sqrt(ty_err_apc.^2 + tx_err_apc.^2);

figure(3);
plot(err_mag_scc, '-ok', 'linewidth', 2);
hold on
plot(err_mag_rpc, '-or', 'linewidth', 2);
plot(err_mag_apc, '-ob', 'linewidth', 2);
hold off
h = legend('SCC', 'RPC', 'APC'); 
set(h, 'fontsize', 16);














