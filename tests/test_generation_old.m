fSize = 12;

region_height = 64;
region_width = 64;

% Center pixels
xc = (region_width  + 1) / 2 + 0.5 * (1 - mod(region_width,  2));
yc = (region_height + 1) / 2 + 0.5 * (1 - mod(region_height, 2));

xv = 1 : region_width;
yv = 1 : region_height;

[XV, YV] = meshgrid(xv - xc, yv - yc);

sx_mean = 3;
sy_mean = 0;

sx_rand = 0;
sy_rand = 0;

sx_min =  3;
sx_max =  3;
sy_min =  0;
sy_max =  0;


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
num_particles = 100;

% %

x_buffer = -10;
y_buffer = -10;

noise_std = 0.0001;

x_min = 1 + x_buffer;
x_max = region_width  - x_buffer;
y_min = 1 + y_buffer;
y_max = region_height - y_buffer;

num_rand_ens = 100;
num_corr_ens = 100;

% Window
g_win = gaussianWindowFilter([region_height, region_width], [0.5, 0.5], 'fraction');
% g_win = rect_lowpass_2D([region_height, region_width], [0.5, 0.5], 'fraction');

f_win = fftshift(fft2(fftshift(g_win)));

f_win_norm = f_win ./ max(f_win(:));

% No window
% g = ones(region_height, region_width);


% Allocate ensemble correlation
cc = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);
pc = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);
cc_abs = zeros(region_height, region_width);



% % Corresponding correlations
for k = 1 : num_corr_ens
    
    fprintf(1, 'Corresponding ensemble %d of %d\n', k, num_corr_ens);
    
    dp_01 = d_mean + d_std * randn(num_particles, 1);
    
    int_01 = d_mean ./ dp_01;
    
    sx_cur = sx_min + (sx_max - sx_min) * rand;
    sy_cur = sy_min + (sy_max - sy_min) * rand;
    
    x_01 = ((x_max - x_min) * rand(num_particles, 1) + x_min);
    y_01 = ((y_max - y_min) * rand(num_particles, 1) + y_min);
    
    x_02 = x_01 + sx_cur + sx_rand * randn(num_particles, 1);
    y_02 = y_01 + sy_cur + sy_rand * randn(num_particles, 1);
    
    noise_mat_01 = noise_std * randn(region_height, region_width);
    noise_mat_02 = noise_std * randn(region_height, region_width);
    
    image_01 = generateParticleImage(region_height, region_width,...
    x_01, y_01, dp_01, int_01) + noise_mat_01;

    image_02 = generateParticleImage(region_height, region_width,...
    x_02, y_02, dp_01, int_01) + noise_mat_02;

    f1 = fftn(image_01 .* g_win, [region_height, region_width]);
    f2 = fftn(image_02 .* g_win, [region_height, region_width]);
    
    cc_current = fftshift(f1 .* conj(f2));
    
    pc_current = phaseOnlyFilter(cc_current);
    
    cc_abs_current = abs(cc_current);
    
    cc = cc + cc_current;
    
    cc_abs = cc_abs + cc_abs_current;
    
    pc = pc + pc_current;
    
%     surf(cc_abs);
%     axis square;
%     xlim([1, region_width]);
%     ylim([1, region_height]);
%     pause(0.1);
    
    
end


% Allocate ensemble correlation
cc_rand = zeros(region_height, region_width) + 1i * zeros(region_height, region_width);

% Random correlations
for k = 1 : num_rand_ens
    
    fprintf(1, 'Random ensemble %d of %d\n', k, num_rand_ens);
    
    dp_rand_01 = d_mean + d_std * randn(num_particles, 1);
    dp_rand_02 = d_mean + d_std * randn(num_particles, 1);
    
    int_rand_01 = d_mean ./ dp_rand_01;
    int_rand_02 = d_mean ./ dp_rand_02;
    
    x_rand_01 = ((x_max - x_min) * rand(num_particles, 1) + x_min);
    y_rand_01 = ((y_max - y_min) * rand(num_particles, 1) + y_min);
    
    x_rand_02 = ((x_max - x_min) * rand(num_particles, 1) + x_min);
    y_rand_02 = ((y_max - y_min) * rand(num_particles, 1) + y_min);
    
    noise_mat_01 = noise_std * randn(region_height, region_width);
    noise_mat_02 = noise_std * randn(region_height, region_width);
    
    
    image_rand_01 = generateParticleImage(region_height, region_width,...
    x_rand_01, y_rand_01, dp_rand_01, int_rand_01) + noise_mat_01;

    image_rand_02 = generateParticleImage(region_height, region_width,...
    x_rand_02, y_rand_02, dp_rand_02, int_rand_02) + noise_mat_02;

    f1_rand = fftn(image_rand_01 .* g_win, [region_height, region_width]);
    f2_rand = fftn(image_rand_02 .* g_win, [region_height, region_width]);
    
    cc_rand = cc_rand + fftshift(f1_rand .* conj(f2_rand));
    
%     surf(abs(cc_rand ./ max(cc_rand(:))), 'linewidth', 0.1);
%     axis off
%     caxis([0, 0.1]);
% %     pause(0.1);
%     title(sprintf('NCC, %d pairs', k), 'fontsize', 12);
%     print(1, '-dpng', '-r200', sprintf('~/Desktop/plots/image_%04d.png', k));
%     
%    
    
    
end

% Noise FFT
noise_ft = fftshift(fft2(noise_mat_01));
noise_ft_mag = abs(noise_ft);

% Filter
COEFF = find_sub_fit(cc_abs, abs(cc_rand));

% Min sub
cc_sub = cc_abs - 1 * abs(cc_rand);

% Gaussian fit
[AMPLITUDE, STD_DEV_Y, STD_DEV_X, YC, XC, ARRAY] = fit_gaussian_2D((abs(cc_sub)));

f_filt_01 = (ARRAY);
f_filt_02 = f_filt_01 - min(f_filt_01(:));
f_filt_norm = f_filt_02 ./ max(f_filt_02(:));



apc_spect = pc .* f_filt_norm;
apc_spatial = fftshift(abs(ifft2(fftshift(apc_spect))));

scc_spatial = fftshift(abs(ifft2(fftshift(cc))));
% 
% subplot(2, 2, 1);
% imagesc(image_01.* g_win);
% axis image
% axis off
% title(sprintf('$\\textrm{Image} \\,, d_p = %0.2f$', d_mean), 'FontSIze', fSize, 'interpreter', 'latex');
% % caxis([0, 2])
% 
% subplot(2, 2, 2);
% 
% cc_rand_norm = abs(cc_rand) ./ max(abs(cc_rand(:)));
% 
% surf((cc_rand_norm), 'linewidth', 0.1);
% axis square;
% xlim([1, region_width]);
% ylim([1, region_height]);
% % zlim([0, 1]);
% % caxis([0, 1]);
% title('NCC', 'FontSize', fSize);
% 
% cc_abs_norm = cc_abs ./ max(cc_abs(:));
% 
% subplot(2, 2, 3);
% surf((cc_abs_norm), 'linewidth', 0.1);
% axis square;
% xlim([1, region_width]);
% ylim([1, region_height]);
% % zlim([0, 1]);
% % caxis([0, 1]);
% title('Ensemble CCC magnitude', 'FontSize', fSize);
% 
% 
% ncc = abs(cc_rand);
% 
% cc_sub_norm = abs(cc_sub) ./ max(cc_sub(:));
% 
% subplot(2, 2, 4);
% surf((cc_sub_norm), 'linewidth', 0.1);
% axis square;
% xlim([1, region_width]);
% ylim([1, region_height]);
% % zlim([0, 1]);
% % caxis([0, 1]);
% title('CCC - NCC', 'FontSize', fSize);

% figure(2);
% subplot(1, 5, 1);
% imagesc(image_01 .* g_win);
% axis image
% 
% pc_real = real(pc);
% subplot(1, 5, 2);
% imagesc(pc_real ./ max(pc_real(:)));
% caxis([-1, 1]);
% axis image;
% 
% 
% subplot(1, 5, 3);
% plot(real(pc(xc, :)), 'linewidth', 2);
% axis square;
% xlim([1, region_width]);
% 
% 
% subplot(1, 5, 4);
% imagesc(angle(pc));
% caxis(pi * [-1, 1]);
% axis image
% 
% subplot(1, 5, 5);
% plot(angle(pc(xc, :)), 'linewidth', 2);
% xlim([1, region_width]);
% axis square


 
 
 
 
% % 

subplot(2, 3, 1)
imagesc((image_01 .* g_win)); 
axis image;
title(sprintf('$\\textrm{Image} \\, , d_p = %0.2f$', d_mean), 'interpreter', 'latex', 'FontSize', fSize);


subplot(2, 3, 2);
surf(abs(cc) ./ max(abs(cc(:))));
set(gca, 'view', [0, 0]);
axis square
xlim([1 region_width]);
ylim([1, region_height]);
zlim([0, 1]);
title('Full complex correlation', 'FontSize', fSize);


% 
subplot(2, 3, 3);
surf(abs(cc_rand ./ max(cc_rand(:))));
set(gca, 'view', [0, 0]);
xlim([1 region_width]);
ylim([1, region_height]);
zlim([0, 1]);
axis square
title('Converged NCC', 'FontSIze', fSize);

subplot(2, 3, 4);
surf(abs(cc_sub) ./ max(abs(cc_sub(:))));
set(gca, 'view', [0, 0]);
% axis image;
xlim([1 region_width]);
ylim([1, region_height]);
title('CC - NCC', 'FontSize', fSize);
axis square

subplot(2, 3, 5);
surf(f_filt_norm);
set(gca, 'view', [0, 0]);
% axis image;
xlim([1 region_width]);
ylim([1, region_height]);
zlim([0, 1]);
title('Gaussian fit', 'FontSize', fSize);
axis square

subplot(2, 3, 6);
surf(apc_spatial ./ max(apc_spatial(:)));
xlim([1 region_width]);
ylim([1, region_height]);
set(gca, 'view', [0, 1]);
% zlim(0.1 * [-1, 1]);
axis square
title('APC plane', 'FontSize', fSize);

% print(1, '-dpng', '-r300', sprintf('~/Desktop/plots_02/plot_dp_%0.2f.png', d_mean));

% 
% subplot(2, 3, 5);
% surf(scc_spatial);
% set(gca, 'view', [0, 1]);
% % axis image;
% xlim([1 region_width]);
% ylim([1, region_height]);
% axis square
% title('SCC');
% 

% figure(2); surf(real(cc_exact));


% figure(2);
% p1 = im_norm_01(y1, :); 
% p2 = im_02_shift(y1, :);
% pd = im_diff(y1, :);
% plot(p1, '-k');
% hold on
% plot(p2, '-r');
% plot(pd, '--k');
% hold off
% axis square
% xlim(x1 + [-10, 10]);
% ylim([-2, 2]);

% dx_std = std(x2 - x1);
% dy_std = std(y2 - y1);
% std_mag = sqrt(dx_std.^2 + dy_std.^2);
% fprintf(1, 'Std = %0.4f\n', std_mag);
% 






