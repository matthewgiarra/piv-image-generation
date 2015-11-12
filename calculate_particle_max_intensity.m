function [particle_max_intensity, particle_image_diameter, emission_power]  = ...
    calculate_particle_max_intensity(MAGNIFICATION, ...
    PARTICLE_DIAMETER, WAVELENGTH_MICRONS, NA, FOCAL_LENGTH, Z);
% Calculate particle image diameter in diffraction limited optics.
% Microscope model taken from Olsen & Adrian 2000

% Airy function constant
B2 = 3.67;

% Objective Magnification
M = MAGNIFICATION;

% Particle diameter
dp = PARTICLE_DIAMETER;

% Focal length
f = FOCAL_LENGTH;

% Wavelength
L = WAVELENGTH_MICRONS;

% F number
f_num = 2 * NA;

% Lens aperture diameter
Da = f / f_num;

% Point response diameter
ds = 2.44 * (M + 1) * L * f_num;

% In focus diameter
de_focused = sqrt(M^2 * dp.^2 + ds^2);

% Emission power constant
emission_power = mean(4 * pi * de_focused.^2 * f^2 ./ (Da^2 * B2));

% Particle image diameter
particle_image_diameter = ...
    calculate_particle_image_diameter(M, dp, L, NA, f, Z);
  
% Max intensity
particle_max_intensity = emission_power * Da^2 * B2 ./...
    (4 * pi * particle_image_diameter.^2 ...
    .* (f + Z).^2);
    
end

