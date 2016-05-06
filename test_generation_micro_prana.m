
addpath ~/Desktop/spectral-phase-correlation/scripts
addpath ~/Desktop/spectral-phase-correlation/jobfiles
addpath ~/Desktop/spectral-phase-correlation/filtering
addpath ~/Desktop/piv-image-generation
addpath ~/Desktop/piv-image-generation/jobfiles

% Region sizes
region_height = 128;
region_width = 128;

% Number of images
num_images = 1;

% Displacements
sx = 3.00 + 1 * 0.36436;
sy = 0;

% Diffusion
diffusion_stdev = 0;

% Mean particle diameter                               
d_mean = 1;

% Particle concentration
c = 2E-2;

% Window
g = gaussianWindowFilter([region_height, region_width], [0.4, 0.4], 'fraction');

% Run the SPC job list
[I1, I2] = microPIVimg(sx, sy);

% Image size
[height, width] = size(I1);

xc = round(width/2);
yc = round(height/2);

xRange = xc : (xc + region_width-1);
yRange = yc : (yc + region_height - 1);

image_01 = I1(yRange, xRange);
image_02 = I2(yRange, xRange);

% FFTs
f1 = fft2(image_01 .* g);
f2 = fft2(image_02 .* g);

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

surf(rpc ./ max(rpc(:)));
title('RPC plane of single images');

% % Do the analysis
% job_save_path_list = runMonteCarloCorrelationJobFile(spc_joblist);
% 
% % Load the results
% load(job_save_path_list{1}{1});

% % Calculate the errors
% tx_err_scc = (TX_TRUE - tx_scc);
% ty_err_scc = (-TY_TRUE - ty_scc);
% err_mag_scc = sqrt(ty_err_scc.^2 + tx_err_scc.^2);
% 
% tx_err_rpc = (TX_TRUE - tx_rpc);
% ty_err_rpc = (-TY_TRUE - ty_rpc);
% err_mag_rpc = sqrt(ty_err_rpc.^2 + tx_err_rpc.^2);
% 
% tx_err_apc = (TX_TRUE - tx_apc);
% ty_err_apc = (-TY_TRUE - ty_apc);
% err_mag_apc = sqrt(ty_err_apc.^2 + tx_err_apc.^2);


% close all;

figure(1);
subplot(1, 2, 1)
imagesc(image_01); axis image
axis off
title(sprintf('$t_x = %0.2f, t_y = %0.2f$', sx, sy), ...
    'interpreter', 'latex', ...
    'FontSize', 16);

subplot(1, 2, 2);
imagesc(fftshift(angle(pc)));
axis image
axis off
title({'Phase Correlation', sprintf('Diffusion = %0.2f', diffusion_stdev)}, ...
    'interpreter', 'latex', ...
    'FontSize', 16);

% p = get(gca, 'position');
% p(1) = 0.52;
% set(gca, 'position', p);
% 
% 
% correlation_plot_save_name = sprintf('ensemble_correlation_plot_h_%d_w_%d_tx_%0.2f_ty_%0.2f_diff_%0.2f.png', ...
%     region_height, region_width, sx, sy, diffusion_stdev);

% correlation_plot_save_dir = '~/Desktop/ensemble_plots/correlation_plots';
% 
% if ~exist(correlation_plot_save_dir, 'dir')
%     mkdir(correlation_plot_save_dir);
% end

% correlation_plot_save_path = fullfile(correlation_plot_save_dir, correlation_plot_save_name);
% print(1, '-dpng', '-r300', correlation_plot_save_path);


% figure(2);
% plot(err_mag_scc, '-k', 'linewidth', 2);
% hold on
% plot(err_mag_rpc, '-r', 'linewidth', 2);
% plot(err_mag_apc, '-b', 'linewidth', 2);
% hold off
% h = legend('SCC', 'RPC', 'APC'); 
% set(h, 'fontsize', 16);
% set(gca, 'FontSize', 16);
% axis square
% xlabel('Number of pairs', 'FontSize', 16);
% ylabel('Translation error magnitude (pixels)', 'FontSize', 16);
% title({['$t_x = ' num2str(sx, '%0.2f') ', t_y = ' ...
%     num2str(sy, '%0.2f') '\, \textrm{pix}$'], ['$ \left(u^{\prime}, v^{\prime}\right) = ' ...
%     num2str(diffusion_stdev, '%0.2f') ' \, \textrm{pix/frame} $'], ...
%     ['$ \bar{d_p} = ' ...
%     num2str(d_mean, '%0.2f') ' \, \textrm{pix} $']
%     } , 'interpreter', 'latex');
% 
% ylim([0, 1]);
% 
% error_plot_save_name = sprintf('ensemble_plot_h_%d_w_%d_tx_%0.2f_ty_%0.2f_diff_%0.2f.png', ...
%     region_height, region_width, sx, sy, diffusion_stdev); 
% 
% error_plot_save_dir = '~/Desktop/ensemble_plots/error_plots';
% 
% if ~exist(error_plot_save_dir, 'dir')
%     mkdir(error_plot_save_dir);
% end
% 
% error_plot_save_path = fullfile(error_plot_save_dir, error_plot_save_name);
% 
% print(2, '-dpng', '-r300', error_plot_save_path);
% 
% 












