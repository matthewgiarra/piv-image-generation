fSize = 12;

region_height = 64;
region_width = 64;

% Image dimensions
image_width  = 1280;
image_height = 1024;

window_fraction = 0.5 * [1, 1];

% Center pixels
xc = (region_width  + 1) / 2 + 0.5 * (1 - mod(region_width,  2));
yc = (region_height + 1) / 2 + 0.5 * (1 - mod(region_height, 2));

xv = 1 : region_width;
yv = 1 : region_height;

[XV, YV] = meshgrid(xv - xc, yv - yc);

% Mean displacements
sx_mean = 3;
sy_mean = 0;

% Random displacements
s_rand = 1;
sx_rand = s_rand;
sy_rand = s_rand;

% Particle positions buffer
x_buffer = -100;
y_buffer = -100;

% Number of non-corresponding correlations for the NCC
num_ens_neq = 50000;

% % Particle stuff

% Standard deviation of particle image diameters
d_std = 0;

% Mean particle diameter
% d_mean = 5 * sqrt(8);
% % % about sqrt(8) pix seems to give the highest SNR? is this some
% optimization?
d_mean = 1 * sqrt(8) ;
% d_mean = 1 ;

% Particle concentration in particles per pixel
particle_concentration = 5E-3;

% Image noise
noise_std = 1E-4;


% % % % % % % % %

% Image coordinates
xc_img = image_width/2;
yc_img = image_width/2;

% Window
g_win = gaussianWindowFilter([region_height, region_width], window_fraction, 'fraction');
% g_win = rect_lowpass_2D([region_height, region_width], [0.5, 0.5], 'fraction');

% % Image parameters

% Compute the total number of particles
num_particles = round(particle_concentration * image_height * image_width);

% Particle diameters
dp       = d_mean + d_std * randn(num_particles, 1);

% Particle intensities for the correlated images
particle_intensities     = d_mean ./ dp;

% Particle position max and min
x_min = 1 + x_buffer;
x_max = image_width  - x_buffer;
y_min = 1 + y_buffer;
y_max = image_height - y_buffer;

% Particle positions (image 1)
x_01 = (x_max - x_min) * rand(num_particles, 1) + x_min;
y_01 = (y_max - y_min) * rand(num_particles, 1) + y_min;

% Velocity function
U_VECT = lambOseenVortexRingVelocityFunction(0, [x_01(:); y_01(:)]);

% Parse the output of the velocity function
% u = U_VECT(1 : num_particles);
% v = U_VECT(num_particles + 1 : end);

u = sx_mean * ones(num_particles, 1);
v = sy_mean * ones(num_particles, 1);

% Particle positions (image 2)
x_02 = x_01 + u + sx_rand * randn(num_particles, 1);
y_02 = y_01 + v + sy_rand * randn(num_particles, 1);

% Noise matrices
noise_mat_01      = noise_std * randn(image_height, image_width);
noise_mat_02      = noise_std * randn(image_height, image_width);

% % Generate the images
%
% % Generate the first image
image_01 = generateParticleImage(image_height, image_width,...
x_01, y_01, dp, particle_intensities) + noise_mat_01;
%
% Generate the second image
image_02 = generateParticleImage(image_height, image_width,...
x_02, y_02, dp, particle_intensities) + noise_mat_02;
% 
%
%

% Grid spacing
grid_spacing_y = region_height;
grid_spacing_x = region_width;
grid_buffer_x = region_width * [1, 1];
grid_buffer_y = region_height * [1, 1];

% Grid the image
[gx, gy] = gridImage([image_height, image_width], ...
    [grid_spacing_y, grid_spacing_x], ...
    grid_buffer_y, grid_buffer_x);

% Extract the subregions (image 1)
region_matrix_01 = extractSubRegions(image_01, ...
    [region_width, region_height], gx, gy);

% Extract the subregions (image 2)
region_matrix_02 = extractSubRegions(image_02, ...
    [region_width, region_height], gx, gy);

% Number of regions
num_regions = size(region_matrix_01, 3);

% Augment the list by a lot so we can get rid of the equal-region pairs
% and still have plenty left over
% List of the first random correlation regions
region_numbers_neq_aug_01 = round((num_regions - 1) * rand(2 * num_ens_neq, 1) + 1);

% List of the second random correlation regions
region_numbers_neq_aug_02 = round((num_regions - 1) * rand(2 * num_ens_neq, 1) + 1);

% Find the equivalent indices
equiv_inds = find(region_numbers_neq_aug_01 == region_numbers_neq_aug_02);

% Pop the equivalent indices
region_numbers_neq_aug_01(equiv_inds) = [];
region_numbers_neq_aug_02(equiv_inds) = [];

% Take the first num_ens_neq of regions
region_numbers_neq_01 = region_numbers_neq_aug_01(1 : num_ens_neq);
region_numbers_neq_02 = region_numbers_neq_aug_02(1 : num_ens_neq);




cc_neq_full = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);
cc_neq_abs_full = zeros(region_height, region_width);
cc_neq_dc = zeros(num_ens_neq, 1);
cc_neq_real_sq = zeros(region_height, region_width);
cc_neq_imag_sq = zeros(region_height, region_width);


% Non-corresponding correlations
for k = 1 : num_ens_neq
    
    % Inform the user
    if mod(k, 10) == 0
        fprintf(1, 'On non-corresponding region %d of %d\n', k, num_ens_neq);
    end
    
    % Pick the region numbers
    region_number_01 = region_numbers_neq_01(k);
    region_number_02 = region_numbers_neq_02(k);
    
    % Extract the region
    region_01 = region_matrix_01(:, :, region_number_01);
    region_02 = region_matrix_02(:, :, region_number_02);
    
    % Fourier transforms of correlated images
    f1 = fftn(region_01 .* g_win, [region_height, region_width]);
    f2 = fftn(region_02 .* g_win, [region_height, region_width]);
    
    % Complex cross correlation of the correlated images
    cc_neq_cur = fftshift(f1 .* conj(f2));
    
    % Running full coresponding correlation
    cc_neq_full = cc_neq_full + cc_neq_cur;
    
    % Running magnitude of the current corresponding CC
    cc_neq_abs_full = cc_neq_abs_full + abs(cc_neq_cur);
  
    % Running DC component of corresponding CC
    cc_neq_dc(k) = real(cc_neq_cur(yc, xc));
    
    % Real and imaginary squared parts
    cc_neq_real_sq = cc_neq_real_sq + (real(cc_neq_cur)).^2;
    cc_neq_imag_sq = cc_neq_imag_sq + (imag(cc_neq_cur)).^2;
    
    % Grid points
    gx_01 = gx(region_number_01);
    gy_01 = gy(region_number_01);
    
    gx_02 = gx(region_number_02);
    gy_02 = gy(region_number_02);
    
%     % Show the image
%     figure(1);
%     imagesc(image_01); 
%     axis image;
%     colormap gray;
%     hold on
%     plot(gx_01, gy_01, 'oy', 'markerfacecolor', 'y', 'markersize', 10);
%     plot(gx_02, gy_02, 'or', 'markerfacecolor', 'r', 'markersize', 10);
%     hold off
%     
%     pause(0.1);

end

% Scaled non-corresponding correlation
cc_neq_scaled = cc_neq_full ./ num_ens_neq;



% Running ensemble DC components
cc_eq_dc = zeros(num_regions, 1);

% Allocate ensemble correlation
cc_eq_full = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);
cc_eq_abs_full  = zeros(region_height, region_width);


cc_eq_dc = zeros(num_regions, 1);
cc_eq_real_sq = zeros(region_height, region_width);
cc_eq_imag_sq = zeros(region_height, region_width);

cc_eq_sub_abs_full = zeros(region_height, region_width);
cc_eq_sub_full = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

% Corresponding-region correlations
for k = 1 : num_regions
    
    % Inform the user
    if mod(k, 10) == 0
        fprintf(1, 'On corresponding region %d of %d\n', k, num_regions);
    end
    
    % Extract the region
    region_01 = region_matrix_01(:, :, k);
    region_02 = region_matrix_02(:, :, k);
    
    % Fourier transforms of correlated images
    f1 = fftn(region_01 .* g_win, [region_height, region_width]);
    f2 = fftn(region_02 .* g_win, [region_height, region_width]);
    
    % Complex cross correlation of the correlated images
    cc_eq_cur = fftshift(f1 .* conj(f2));
    
    % Effective number of particles
    np_eff = sqrt(abs(real(cc_eq_cur(yc, xc))));
    
    % Scaling factor 
    
    % Subtracted
    cc_eq_sub = cc_eq_cur - np_eff * (np_eff - 1) * cc_neq_scaled ./ max(real(cc_neq_scaled(:)));
    
    % Full subtracted
    cc_eq_sub_full = cc_eq_sub_full + cc_eq_sub;
    
    % Running full coresponding correlation
    cc_eq_full = cc_eq_full + cc_eq_cur;
    
    % Running magnitude of the current corresponding CC
    cc_eq_abs_full = cc_eq_abs_full + abs(cc_eq_cur);
  
    % Running DC component of corresponding CC
    cc_eq_dc(k)  = real(cc_eq_cur(yc, xc));
    
    % Real and imaginary squared parts
    cc_eq_real_sq = cc_eq_real_sq + (real(cc_eq_cur)).^2;
    cc_eq_imag_sq = cc_eq_imag_sq + (imag(cc_eq_cur)).^2;
    
    cc_eq_sub_abs_full = cc_eq_sub_abs_full + abs(cc_eq_sub);

end

surf(cc_eq_sub_abs_full); 
axis square;
set(gca, 'view', [0, 0]);


% 
% 
% cc_eq_mag = sqrt(cc_eq_real_sq + cc_eq_imag_sq);
% cc_neq_mag = sqrt(cc_neq_real_sq + cc_neq_imag_sq);
% 
% cc_mag_sub = cc_eq_mag - (num_regions / num_ens_neq * cc_neq_mag);
% 
% 
% cc_sub = cc_eq_full - (num_regions / num_ens_neq * cc_neq_full);
% 
% cc_sub_abs = abs(cc_sub);
% 
% cc_abs_sub = cc_eq_abs_full - (num_regions / num_ens_neq * abs(cc_neq_full));
% 
% cc_spatial_raw = fftshift(abs(ifft2(fftshift(cc_eq_full))));
% 
% cc_sub_spatial = fftshift(abs(ifft2(fftshift(cc_sub))));
% 

% subplot(2, 2, 1);
% surf(real(cc_eq_full));
% set(gca, 'view', [0, 0.1]);
% axis square;
% xlim([1, region_width]);
% ylim([1, region_height]);
% 
% subplot(2, 2, 2);
% surf(real(cc_neq_full));
% axis square;
% xlim([1, region_width]);
% ylim([1, region_height]);
% 
% subplot(2, 2, 3);
% % imagesc(real(cc_sub) ./ max(abs(real(cc_sub(:)))));
% surf(real(cc_sub) ./ max(abs(real(cc_sub(:)))));
% % caxis([-1, 1]);
% set(gca, 'view', [0, 0]);
% axis square;
% xlim([1, region_width]);
% ylim([1, region_height]);
% 
% 
% subplot(2, 2, 4);
% surf(abs(cc_abs_sub));
% axis square;
% set(gca, 'view', [0, 0]);
% xlim([1, region_width]);
% ylim([1, region_height]);

% [AMPLITUDE, STD_DEV_Y, STD_DEV_X, YC, XC, ARRAY] = fit_gaussian_2D(abs(cc_sub_abs));
% 
% fprintf('Std X = %0.2f\tStd Y = %0.2f\n', STD_DEV_X, STD_DEV_Y);

% fprintf('Eq DC: %0.2f\t NEQ DC: %0.2f\n', cc_eq_dc, cc_neq_dc);

% n_eq =

% 
% abs_running_diff = max(cc_eq_abs_full(:)) - max(cc_neq_abs_full(:));
% 
% n = length(cc_eq_dc_running) - 1;
% dc_running_diff = (cc_eq_dc_running - cc_neq_dc_running)/n;
% 
% % cc_sub = cc_eq_abs_full - pct * abs(cc_neq_full);
% 
% N_eff_eq = sqrt(max(cc_eq_abs_full(:)));
% N_eff_neq = sqrt(max(abs(cc_neq_full(:))));
% 
% pct = N_eff_neq / N_eff_eq;
% 
% 
% 
% % cc_sub = cc_eq_abs_full - pct * abs(cc_neq_full);
% % cc_fit = 
% % cc_sub = cc_eq_abs_full -  abs(cc_neq_full);
% cc_sub = cc_eq_abs_full - cc_neq_abs_full;
% 
% subplot(1, 2, 1);
% plot(cc_eq_dc_running, 'ok');
% hold on
% plot(cc_neq_dc_running, 'or');
% hold off
% axis square
% h = legend('EQ', 'NEQ');
% set(h, 'fontsize', 16);
% 
% subplot(1, 2, 2);
% surf(cc_sub);
% axis square
% set(h, 'fontsize', 16);
% title('Difference', 'FontSize', 16);
% set(gca, 'view', [0, 0]);
% 
% % Real part of ensemble corresponding correlation
% rcc = real(cc_eq_full);
% 
% % Real part of ensemble non-corresponding correlation
% ncc = real(cc_neq_full);

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







