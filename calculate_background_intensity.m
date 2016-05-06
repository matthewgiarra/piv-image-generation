function background_intensity = calculate_background_intensity(...
    CHANNEL_DEPTH_MICRONS, PARTICLE_CONCENTRATION, PARTICLE_DIAMETER_MICRONS, ...
    WAVELENGTH_MICRONS, MAGNIFICATION, NA, WORKING_DISTANCE_MICRONS, ...
    FOCAL_LENGTH_MICRONS, OBJECT_PLANE_DISTANCE_MICRONS)

% This function calculates the intensity of the background illumination
% in a micro PIV image.
% 
% Source: Reduced Equation 9 in Olsen & Adrian 2000

% "Adrian and Yao found that for a Gaussian distribution
% to best approximate an Airy distribution, B^2 = 3.67"
% Olsen and Adrian 2000 equation 5
B2 =  3.67;

% Channel depth (microns)
L = CHANNEL_DEPTH_MICRONS;

% Object plane distance
so = OBJECT_PLANE_DISTANCE_MICRONS;

% Concentration (particles per micron^3)
C = PARTICLE_CONCENTRATION;

% Particle diameter in microns
dp = PARTICLE_DIAMETER_MICRONS;

% Magnification
M = MAGNIFICATION;

% Focal length
f = FOCAL_LENGTH_MICRONS;

% Aperture diameter
Da_microns = 2 * NA * FOCAL_LENGTH_MICRONS;

% Working distance in microns
wd = WORKING_DISTANCE_MICRONS;

% Diffraction limited diameter
ds = 2.44 * (M + 1) * WAVELENGTH_MICRONS / (2 * NA);

% In focus particle diameter
de = sqrt(M^2 * dp^2 + ds^2);

% Light power from a particle
Jp = 4 * pi * de.^2 * WORKING_DISTANCE_MICRONS^2 / (B2 * Da_microns^2);

% Background intensity
background_intensity = Jp * C  * L * Da_microns^2 / (16 * MAGNIFICATION^2 * (wd - so) * (wd - so + L));

end

