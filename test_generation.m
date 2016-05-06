region_height = 128;
region_width = 128;

sx = 5.0 + 1 * rand;
sy = 3.0 + 1 * rand;

num_particles = 50;

x_buffer = -20;
y_buffer = -20;

x_min = 1 + x_buffer;
x_max = region_width - x_buffer;
y_min = 1 + y_buffer;
y_max = region_height - y_buffer;


% Window
g = gaussianWindowFilter([region_height, region_width], [0.3, 0.3], 'fraction');

% No window
% g = ones(region_height, region_width);

% Standard deviation of particle image diameters
d_std = 3;

% Mean particle diameter
% d_mean = 5 * sqrt(8);
d_mean = 1 * sqrt(8) ;

% All the particle diameters
particle_diameters = d_mean + d_std * randn(num_particles, 1);
% particle_diameters = d_mean;

% Particle intensities
particle_intensities = d_mean ./ (particle_diameters);

x1 = round((x_max - x_min) * rand(num_particles, 1) + x_min);
y1 = round((y_max - y_min) * rand(num_particles, 1) + y_min);

[X, Y] = meshgrid(1 : region_width, 1 : region_height);

z = ones(size(x1));

% x1 = 64;
% y1 = 64;
% 
x2 = x1 + sx;
y2 = y1 + sy;

% std_dev = particle_diameters 

% image_01 = exp(-(X - x1).^2 / (2 * d_mean^2)) .* exp(-(Y - y1).^2 / (2 * d_mean^2));
% image_02 = exp(-(X - x2).^2 / (2 * d_mean^2)) .* exp(-(Y - y2).^2 / (2 * d_mean^2));
% image_02 = circshift(image_01, [sy, sx]);

image_01 = generateParticleImage(region_height, region_width,...
    x1, y1, particle_diameters, particle_intensities);



image_02 = generateParticleImage(region_height, region_width,...
    x2, y2, particle_diameters, particle_intensities);


% image_01_micro = ...
% generate_micro_piv_image(x1, y1, z, ...
%     region_height, region_width, d_mean ,...
%     pixel_size_microns, particle_concentration, channel_depth_microns, ...
%     objective_magnification, focal_length_microns, NA, ...
%     wavelength_microns, intensity_fraction);



im_norm_01 = image_01 ./ max(image_01(:));
im_norm_02 = image_02 ./ max(image_02(:));

% RPC Filter
rpc_filter = spectralEnergyFilter(region_height, region_width, d_mean);

% % Pixel difference
% im_02_shift = circshift(im_norm_02, [-sy, -sx]);
% im_diff = (im_02_shift - im_norm_01);


imagesc(image_01); axis image; 

f1 = fftn(image_01 .* g, [region_height, region_width]);
f2 = fftn(image_02 .* g, [region_height, region_width]);

cc = f1 .* conj(f2);

pc = phaseOnlyFilter(cc);

rpc_spect = pc .* fftshift(rpc_filter);

scc = fftshift(abs(real(ifft2(cc))));
gcc = fftshift(abs(real(ifft2(pc))));
rpc = fftshift(abs(real(ifft2(rpc_spect))));


subplot(2, 2, 1)
imagesc((image_01)); 
axis image;

% subplot(1, 3, 2);
% imagesc(im_diff); axis image;
% axis off

subplot(2, 2, 2);
imagesc(fftshift(angle(pc)));
axis image
axis off

subplot(2, 2, 3);
mesh(scc ./ max(scc(:)), 'edgecolor', 'black');

subplot(2, 2, 4);
mesh(rpc ./ max(rpc(:)), 'edgecolor', 'black');

scc_plane = SCC(image_01 .* g, image_02.*g, 1);


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










