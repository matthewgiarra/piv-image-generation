%#codegen
function IMAGEOUT = generateParticleImage(HEIGHT, WIDTH, X, Y, ...
    PARTICLE_DIAMETERS, PARTICLE_MAX_INTENSITIES);
% This function generate a synthetic PIV particle image.
% Usage: IMAGEOUT = generateParticleImage(HEIGHT, WIDTH, X, Y, PARTICLEDIAMETERS, PARTICLEMAXINTENSITIES)
%
% INPUTS:
%   HEIGHT = Number of rows in the image to generate (integer).
%   WIDTH  = Number of columns in the image to generate (integer).
%   X = Column vector containing the non-integer column
%       positions of the particles to generate.
%   Y = Column vector containing the non-integer row
%       positions of the particles to generate.
%   PARTICLEDIAMETERS = Column vector containing the non-integer
%       diameters of the particles to generate.
%   PARTICLEMAXINTENSITIES = Column vector containing the maximum
%       intensities of all the particles to generate.
%
% OUTPUTS:
%   IMAGEOUT = [HEIGHT x WIDTH] matrix containing the rendered particle
%   image.

% Create a placeholder to store the generated image
imagePlaceholder = zeros(HEIGHT, WIDTH); 

% Only render particles whose intensities are above some threshold,
% which we specify here as the value of a zero-mean Gaussian function
% at two standard deviations (so that the "width" of a symmetric 
% Gaussian beam is four times the beam's standard deviation).
cutoff_intensity = exp(-2);

% Determine which particles to render. Only render particles whose images lie
% within the image domain.
particlesToRender = X <= WIDTH & X >= 1 & Y <= HEIGHT & Y >= 1 & ...
    PARTICLE_MAX_INTENSITIES >= cutoff_intensity;

% Extract the coordinates and brightness of the particles to be rendered.
xRender = X(particlesToRender);
yRender = Y(particlesToRender);
particleIntensities = PARTICLE_MAX_INTENSITIES(particlesToRender);
particleDiameters = PARTICLE_DIAMETERS(particlesToRender);

% Determine the number of particles to render
numberOfParticlesToRender = length(xRender);

% Determine the miniumum and maximum columns (leftmost and rightmost pixels) in the image
% to which each particle contributes some intensity,
% fractional values
minRenderedCols_fract = max(1,     (xRender - 0.5 * particleDiameters));
maxRenderedCols_fract = min(WIDTH, (xRender + 0.5 * particleDiameters));

% Determine the minimum and maximum rows (topmost and bottommost pixels) in
% the image to which each particle contributes some intensity, 
% fractional values
minRenderedRows_fract = max(1,      (yRender - 0.5 * particleDiameters));
maxRenderedRows_fract = min(HEIGHT, (yRender + 0.5 * particleDiameters));

% Take the whole-number part of the fractional values.
minRenderedCols = minRenderedCols_fract - rem(minRenderedCols_fract, 1);
maxRenderedCols = maxRenderedCols_fract - rem(maxRenderedCols_fract, 1);

% Take the whole-number part of the fractional values.
minRenderedRows = minRenderedRows_fract - rem(minRenderedRows_fract, 1);
maxRenderedRows = maxRenderedRows_fract - rem(maxRenderedRows_fract, 1);

% Square root of 8; just calculate this once
sqrt8 = sqrt(8);

% Generate the intensities from each particle.
for p = 1 : numberOfParticlesToRender
    for c = minRenderedCols(p) : maxRenderedCols(p)
        for r = minRenderedRows(p) : maxRenderedRows(p)
           imagePlaceholder(r, c) = imagePlaceholder(r, c) + ...
               particleIntensities(p) * (particleDiameters(p))^2 * pi / 32 *...
                   (erf( sqrt8 *  (c - xRender(p) + 0.5) ...
                   / particleDiameters(p)) - erf(sqrt8 * ...
                   (c - xRender(p) - 0.5) / particleDiameters(p))) * ...
                   ...
                   (erf( sqrt8 *  (r - yRender(p) + 0.5) ...
                   / particleDiameters(p)) - erf(sqrt8 * ...
                   (r - yRender(p) - 0.5) / particleDiameters(p)));
            
        end
    end    
end

% Keep the image from over-saturating by limiting its maximum value to one.
if max(imagePlaceholder(:)) > 1
    imagePlaceholder = imagePlaceholder ./ max(imagePlaceholder(:));
end

% Flip the image, which contains no noise
IMAGEOUT = flipud(imagePlaceholder);

end
