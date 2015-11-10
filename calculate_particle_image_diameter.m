function particle_image_diameter = calculate_particle_image_diameter(...
   MAGNIFICATION, PARTICLE_DIAMETER, WAVELENGTH_MICRONS, NA, FOCAL_LENGTH, Z);
% Calculate particle image diameter in diffraction limited optics.
% Microscope model taken from Olsen & Adrian 2000

% Magnification
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

% Image diameter
particle_image_diameter = sqrt(M^2 * dp ^2 + 5.95 * (M + 1)^2 * L^2 * f_num^2 + ...
    (M * Z * Da).^2 ./ (f + Z).^2);



end