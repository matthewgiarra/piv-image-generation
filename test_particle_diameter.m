M = 10;
dp = 1.00;
wavelength_microns = 0.532;
NA = 0.28;
focal_length_microns = 2E4;
channel_depth_microns = 50;
particle_concentration = 1E-3;
render_intensity_fraction = 0.5;

% Font size
fSize = 16;

% Max Z position to render
z_max = particle_render_depth(...
    channel_depth_microns, particle_concentration, M,...
    dp, wavelength_microns,...
    NA, render_intensity_fraction);

% Vector of positions
z = (linspace(0, z_max, 1000))';

[particle_max_intensity, particle_image_diameter]  = ...
    calculate_particle_image_size(M, ...
    dp, wavelength_microns, NA, focal_length_microns, z);


% Plot particle image diameter
% Particle image diameter
subplot(1, 2, 1);
plot(z , particle_image_diameter, '-k');
xlabel('Z position (microns)', 'FontSize', fSize);
ylabel('Particle image diameter (microns)', 'FontSize', fSize);
title('Diameter', 'FontSize', fSize);
grid on
axis square;
set(gca, 'FontSize', fSize);

subplot(1, 2, 2);
plot(z , particle_max_intensity, '-k');
xlabel('Z position (microns)', 'FontSize', fSize);
ylabel('Normalized max intensity', 'FontSize', fSize);
title('Intensity', 'FontSize', fSize);
grid on
axis square;
set(gca, 'FontSize', fSize);
set(gcf, 'color', 'white');