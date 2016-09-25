fSize = 12;

region_height = 64;
region_width = 64;

window_fraction = 0.3 * [1, 1];

% Center pixels
xc = (region_width  + 1) / 2 + 0.5 * (1 - mod(region_width,  2));
yc = (region_height + 1) / 2 + 0.5 * (1 - mod(region_height, 2));

xv = 1 : region_width;
yv = 1 : region_height;

[XV, YV] = meshgrid(xv - xc, yv - yc);

sx_rand = 100;
sy_rand = 0;

sx_mean = 3;
sy_mean = 0;

sx_range = sx_mean * [1, 1];
sy_range = sy_mean * [1, 1];

% % Particle stuff

% Standard deviation of particle image diameters
d_std = 0;

% Mean particle diameter
% d_mean = 5 * sqrt(8);
% % % about sqrt(8) pix seems to give the highest SNR? is this some
% optimization?
d_mean = 1 * sqrt(8) ;
% d_mean = 1 ;

% Number of particles
num_particles = 500;

% %

x_buffer = -100;
y_buffer = -100;

noise_std = 1E-4;

num_corr_ens = 5000;



% % % % % % % % %

% Window
g_win = gaussianWindowFilter([region_height, region_width], window_fraction, 'fraction');
% g_win = rect_lowpass_2D([region_height, region_width], [0.5, 0.5], 'fraction');

% Max and min particle locations
x_min = 1 + x_buffer;
x_max = region_width  - x_buffer;
y_min = 1 + y_buffer;
y_max = region_height - y_buffer;

% Fourier transform of window
f_win = fftshift(fft2(fftshift(g_win)));

% Normalized fourier transform of window
f_win_norm = f_win ./ max(f_win(:));

% No window
% g = ones(region_height, region_width);

% Allocate ensemble correlation
cc_eq_full = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);
cc_neq_full = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);
cc_eq_abs_full  = zeros(region_height, region_width);
cc_neq_abs_full = zeros(region_height, region_width);
cc_eq_dc = 0;
cc_neq_dc = 0;

% Displacement ranges
sx_min =  sx_range(1);
sx_max =  sx_range(2);
sy_min =  sy_range(1);
sy_max =  sy_range(2);

% Running ensemble DC components
cc_eq_dc_running = zeros(num_corr_ens, 1);
cc_neq_dc_running = zeros(num_corr_ens, 1);

% % Corresponding correlations
for k = 1 : num_corr_ens
    
    if mod(k, 100) == 0
        fprintf(1, 'Corresponding ensemble %d of %d\n', k, num_corr_ens);
    end
    % Particle diameters for the correlated images
    dp_01       = d_mean + d_std * randn(num_particles, 1);
    
    % Particle diameters for the uncorrelated images
    dp_neq_01 = d_mean + d_std * randn(num_particles, 1);
    dp_neq_02 = d_mean + d_std * randn(num_particles, 1);
    dp_neq_03 = d_mean + d_std * randn(num_particles, 1);
    
    % Particle intensities for the correlated images
    int_01      = d_mean ./ dp_01;
    
    % Particle intensities for the uncorrelated images
    int_neq_01 = d_mean ./ dp_neq_01;
    int_neq_02 = d_mean ./ dp_neq_02;
    
    % Displacements for the correlated images
    sx_cur = sx_min + (sx_max - sx_min) * rand;
    sy_cur = sy_min + (sy_max - sy_min) * rand;
    
    % Particle positions (correlated image 1)
    x_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;
    
    % Particle positions (correlated image 2)
    x_02 = x_01 + sx_cur + sx_rand * randn(num_particles, 1);
    y_02 = y_01 + sy_cur + sy_rand * randn(num_particles, 1);
    
    % Particle positions (uncorrelated image 1)
    x_neq_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_neq_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;
    
    % Particle positions (uncorrelated image 2)
    x_neq_02 = (x_max - x_min) * rand(num_particles, 1) + x_min;
    y_neq_02 = (y_max - y_min) * rand(num_particles, 1) + y_min;
    
    % % Noise
    % 
    % Noise for correlated images
    noise_mat_01      = noise_std * randn(region_height, region_width);
    noise_mat_02      = noise_std * randn(region_height, region_width);
    %
    % Noise for uncorrelated images.
    noise_mat_neq_01 = noise_std * randn(region_height, region_width);
    noise_mat_neq_02 = noise_std * randn(region_height, region_width);
    
    % Generate the first correlated image
    image_01 = generateParticleImage(region_height, region_width,...
    x_01, y_01, dp_01, int_01) + noise_mat_01;

    % Generate the second correlated image
    image_02 = generateParticleImage(region_height, region_width,...
    x_02, y_02, dp_01, int_01) + noise_mat_02;

    % Generate the first uncorrelated image
    image_neq_01 = generateParticleImage(region_height, region_width,...
    x_neq_01, y_neq_01, dp_neq_01, int_neq_01) + noise_mat_neq_01;

    % Generate the second uncorrelated image
    image_neq_02 = generateParticleImage(region_height, region_width,...
    x_neq_02, y_neq_02, dp_neq_02, int_neq_02) + noise_mat_neq_02;
   
    % Fourier transforms of correlated images
    f1 = fftn(image_01 .* g_win, [region_height, region_width]);
    f2 = fftn(image_02 .* g_win, [region_height, region_width]);
    
    % Fourier transforms of the uncorrelated images
    f1_neq = fftn(image_neq_01 .* g_win, [region_height, region_width]);
    f2_neq = fftn(image_neq_02 .* g_win, [region_height, region_width]);
    
    % Complex cross correlation of the correlated images
    cc_eq_cur = fftshift(f1 .* conj(f2));
    
    % Complex cross correlation of the un-correlated images
    cc_neq_cur = fftshift(f1_neq .* conj(f2_neq));
    
    % Running full coresponding correlation
    cc_eq_full = cc_eq_full + cc_eq_cur;
    
    % Running full non-oresponding correlation
    cc_neq_full = cc_neq_full + cc_neq_cur;
    
    % Running magnitude of the current corresponding CC
    cc_eq_abs_full = cc_eq_abs_full + abs(cc_eq_cur);
    
    % Running magnitude of the current non-corresponding CC
    cc_neq_abs_full = cc_neq_abs_full + abs(cc_neq_cur);
  
    % Running DC component of corresponding CC
    cc_eq_dc  = cc_eq_dc + real(cc_eq_cur(yc, xc));
    
    % Running DC component of non-corresponding CC
    cc_neq_dc = cc_neq_dc + real(cc_neq_cur(yc, xc));
    
    % Sum of running DC component for corresponding correlation
    cc_eq_dc_running(k) = sum(cc_eq_dc);
    
     % Sum of running DC component for corresponding correlation
    cc_neq_dc_running(k) = sum(cc_neq_dc);
   
end

% fprintf('Eq DC: %0.2f\t NEQ DC: %0.2f\n', cc_eq_dc, cc_neq_dc);

% n_eq =


abs_running_diff = max(cc_eq_abs_full(:)) - max(cc_neq_abs_full(:));

n = length(cc_eq_dc_running) - 1;
dc_running_diff = (cc_eq_dc_running - cc_neq_dc_running)/n;

% cc_sub = cc_eq_abs_full - pct * abs(cc_neq_full);

N_eff_eq = sqrt(max(cc_eq_abs_full(:)));
N_eff_neq = sqrt(max(abs(cc_neq_full(:))));

pct = N_eff_neq / N_eff_eq;



% cc_sub = cc_eq_abs_full - pct * abs(cc_neq_full);
% cc_fit = 
% cc_sub = cc_eq_abs_full -  abs(cc_neq_full);
cc_sub = cc_eq_abs_full - cc_neq_abs_full;

subplot(1, 2, 1);
plot(cc_eq_dc_running, 'ok');
hold on
plot(cc_neq_dc_running, 'or');
hold off
axis square
h = legend('EQ', 'NEQ');
set(h, 'fontsize', 16);

subplot(1, 2, 2);
surf(cc_sub);
axis square
set(h, 'fontsize', 16);
title('Difference', 'FontSize', 16);
set(gca, 'view', [0, 0]);

% Real part of ensemble corresponding correlation
rcc = real(cc_eq_full);

% Real part of ensemble non-corresponding correlation
ncc = real(cc_neq_full);

% % Gaussian fit
% [AMPLITUDE, STD_DEV_Y, STD_DEV_X, YC, XC, ARRAY] = fit_gaussian_2D((abs(cc_sub)));
% 
% f_filt_01 = (ARRAY);
% f_filt_02 = f_filt_01 - min(f_filt_01(:));
% f_filt_norm = f_filt_02 ./ max(f_filt_02(:));
% 
% 
% apc_spect = pc_full .* f_filt_norm;
% apc_spatial = fftshift(abs(ifft2(fftshift(apc_spect))));
% 
% scc_spatial = fftshift(abs(ifft2(fftshift(ph_eq_full))));


% 
% figure(3); 
% subplot(1, 3, 1);
% surf(rcc./ max(rcc(:)));
% zlim([-0.5, 1]);
% xlim([1, region_width]);
% ylim([1, region_height]);
% set(gca, 'view', [0, 0]);
% axis square
% title('CC', 'fontsize', 20);
% 
% subplot(1, 3, 2);
% surf(ncc ./ max(ncc(:)));
% zlim([-0.5, 1]);
% xlim([1, region_width]);
% ylim([1, region_height]);
% set(gca, 'view', [0, 0]);
% axis square
% title('NCC', 'fontsize', 20);
% 
% subplot(1, 3, 3);
% surf(cc_sub ./ max(cc_sub(:)));
% zlim(1.1 * [-1, 1]);
% xlim([1, region_width]);
% ylim([1, region_height]);
% set(gca, 'view', [0, 0]);
% axis square
% title('CC_sub', 'fontsize', 20);
%  







