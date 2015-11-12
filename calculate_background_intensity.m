function background_intensity = calculate_background_intensity(...
    CHANNEL_DEPTH_MICRONS, PARTICLE_CONCENTRATION, MAGNIFICATION,...
    NA, EMISSION_POWER)

% This function calculates the intensity of the background illumination
% in a micro PIV image.
% 
% Source: Reduced Equation 9 in Olsen & Adrian 2000

% Channel depth (microns)
L = CHANNEL_DEPTH_MICRONS;

% Concentration (particles per micron^3)
C = PARTICLE_CONCENTRATION;

% Magnification
M = MAGNIFICATION;

% Background intensity
background_intensity = ...
   EMISSION_POWER * C * L / (4 * M^2 * NA^2);


end