function z_max = particle_render_depth(...
    CHANNEL_DEPTH_MICRONS, PARTICLE_CONCENTRATION, MAGNIFICATION,...
    PARTICLE_DIAMETER_MICRONS, WAVELENGTH_MICRONS,...
    NA, RENDER_INTENSITY_FRACTION);
% This function calculates the maximum distance from the object plane
% at which particles are rendered.
% This criterion is based on the particle maximum intensity being 
% greater than the average background intensity of the image.


% % Channel depth (microns)
% L = 50;
% 
% % Concentration (particles per micron^3)
% C = 0.01;
% 
% % Magnification
% M = 10;
% 
% % Particle diameter (microns)
% dp = 1.00;
% 
% % Illuminating wavelength (microns)
% lambda = 0.532;
% 
% % Numerical aperture
% NA = 0.25;
% 
% % Fraction above background glow
% % intensity that particles must be
% % to be visible (i.e., rendered as particles)
% g = 0.1;

% Channel depth (microns)
L = CHANNEL_DEPTH_MICRONS;

% Concentration (particles per micron^3)
C = PARTICLE_CONCENTRATION;

% Magnification
M = MAGNIFICATION;

% Particle diameter (microns)
dp = PARTICLE_DIAMETER_MICRONS;

% Illuminating wavelength (microns)
lambda = WAVELENGTH_MICRONS;

% Fraction above background glow
% intensity that particles must be
% to be visible (i.e., rendered as particles)
g = RENDER_INTENSITY_FRACTION;

% f number of lens
f = 2 * NA;

% Constant that makes a Gaussian best approximate
% an Airy distribution (Adrian and Yao, 1985) B^2 = 3.67
B2 = 3.67;

% Max distance from focal plane to render particles
z_max = f * ...
    (4 * B2 / (g * pi * C * L) - ...
    1 / M^2 * (M^2 * dp^2 + ...
    5.95 * (M + 1)^2 * lambda^2 * f^2))^(1/2);


end