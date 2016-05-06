function [particle_max_intensity, particle_image_diameters]  = ...
    calculate_particle_max_intensity(MAGNIFICATION, ...
    PARTICLE_DIAMETER, WAVELENGTH_MICRONS, NA, Z);
% Calculate particle image diameter in diffraction limited optics.
% Microscope model taken from Olsen & Adrian 2000

% Airy function constant
B2 = 3.67;

% Objective Magnification
M = MAGNIFICATION;

% Particle diameter
dp = PARTICLE_DIAMETER;

% Wavelength
L = WAVELENGTH_MICRONS;

% Point response diameter
ds = 2.44 * (M + 1) * L / (2 * NA);

% In focus diameter
de_focused = sqrt(M^2 * dp.^2 + ds^2);

% Particle image diameter
particle_image_diameters = ...
    calculate_particle_image_diameter(M, dp, L, NA, Z);
  
% Max intensity
particle_max_intensity = (de_focused ./ particle_image_diameters).^2;
    
end

