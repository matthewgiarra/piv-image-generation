region_height = 128;
region_width = 128;

sx = 5 + 1 * rand;
sy = 3 + 1 * rand;

num_particles = 500;

x_buffer = -20;
y_buffer = -20;

x_min = 1 + x_buffer;
x_max = region_width - x_buffer;
y_min = 1 + y_buffer;
y_max = region_height - y_buffer;

% Window
g = gaussianWindowFilter([region_height, region_width], [0.5, 0.5], 'fraction');

% No window
% g = ones(region_height, region_width);

% Standard deviation of particle image diameters
d_std = 10;

% Mean particle diameter
d_mean = 2 * sqrt(8);

% All the particle diameters
particle_diameters = d_mean + d_std * randn(num_particles, 1);

% Particle intensities
particle_intensities = d_mean ./ (particle_diameters);

x1 = (x_max - x_min) * rand(num_particles, 1) + x_min;
y1 = (y_max - y_min) * rand(num_particles, 1) + y_min;

x2 = x1 + sx;
y2 = y1 + sy;

image_01 = generateParticleImage(region_height, region_width,...
    x1, y1, particle_diameters, particle_intensities) .* g;

image_02 = generateParticleImage(region_height, region_width,...
    x2, y2, particle_diameters, particle_intensities) .* g;

imagesc(image_01); axis image; 

f1 = fft2(image_01);
f2 = fft2(image_02);

cc = f1 .* conj(f2);

pc = phaseOnlyFilter(cc);

subplot(1, 2, 1)
imagesc(image_01); axis image
axis off

subplot(1, 2, 2);
imagesc(fftshift(angle(pc)));
axis image
axis off













