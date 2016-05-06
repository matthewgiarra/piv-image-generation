function background_intensity = calculate_background_intensity(...
    CHANNEL_DEPTH_MICRONS, PARTICLE_CONCENTRATION, ...
    MAGNIFICATION, NA, OBJECT_PLANE_DISTANCE_MICRONS, FOCAL_LENGTH_MICRONS)

% This function calculates the intensity of the background illumination
% in a micro PIV image.
% 
% Source: Reduced Equation 9 in Olsen & Adrian 2000

% Channel depth (microns)
L = CHANNEL_DEPTH_MICRONS;

% Object plane distance
so = OBJECT_PLANE_DISTANCE_MICRONS;

% Concentration (particles per micron^3)
C = PARTICLE_CONCENTRATION;

% Magnification
M = MAGNIFICATION;

% Focal length
f = FOCAL_LENGTH_MICRONS;

% Background intensity
background_intensity = C * L * f^2 * NA^2 / (M * (4 * so^2 - L^2));


end