function particle_image_diameter_microns = calculate_particle_image_diameter(...
   MAGNIFICATION, PARTICLE_DIAMETER_MICRONS, WAVELENGTH_MICRONS, NA, Z_MICRONS);
% Calculate particle image diameter in diffraction limited optics.
% Microscope model taken from Olsen & Adrian 2000

% Magnification
M = MAGNIFICATION;

% Particle diameter
dp = PARTICLE_DIAMETER_MICRONS;

% Wavelength
L = WAVELENGTH_MICRONS;

% Image diameter
particle_image_diameter_microns =...
    sqrt(M^2 * dp .^2 + 5.95 * (M + 1)^2 * L^2 /  (4 * NA^2) + ...
   (2 * NA * M * Z_MICRONS).^2);

end