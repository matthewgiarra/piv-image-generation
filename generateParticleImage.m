%#codegen
function IMAGEOUT = generateParticleImage(HEIGHT, WIDTH, X, Y, ...
    PARTICLEDIAMETERS, PARTICLEMAXINTENSITIES);
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

% Determine which particles to render. Only render particles whose images lie
% within the image domain. 
particlesToRender = X <= WIDTH & X >= 1 & Y <= HEIGHT & Y >= 1;

% Extract the coordinates and brightness of the particles to be rendered.
xRender = X(particlesToRender);
yRender = Y(particlesToRender);
particleIntensities = PARTICLEMAXINTENSITIES(particlesToRender);
particleDiameters = PARTICLEDIAMETERS(particlesToRender);

% Determine the number of particles to render
numberOfParticlesToRender = length(xRender);

% Average particle diameter
mean_particle_diameter = mean(particleDiameters);

% Determine the miniumum and maximum columns (leftmost and rightmost pixels) in the image
% to which each particle contributes some intensity.
minRenderedCols = max(1,    floor(xRender - 0.75 * mean_particle_diameter));
maxRenderedCols = min(WIDTH, ceil(xRender + 0.75 * mean_particle_diameter));

% Determine the minimum and maximum rows (topmost and bottommost pixels) in
% the image to which each particle contributes some intensity.
minRenderedRows = max(1,     floor(yRender - 0.75 * mean_particle_diameter));
maxRenderedRows = min(HEIGHT, ceil(yRender + 0.75 * mean_particle_diameter));

% Generate the intensities from each particle.
for p = 1 : numberOfParticlesToRender
    for c = minRenderedCols(p) : maxRenderedCols(p)
        for r = minRenderedRows(p) : maxRenderedRows(p)
            imagePlaceholder(r, c) = imagePlaceholder(r, c) + ...
                particleIntensities(p) * (particleDiameters(p))^2 * pi / 32 *...
                   (erf( particleDiameters(p) *  (c - xRender(p) + 0.5) / particleDiameters(p)) - erf(particleDiameters(p) * (c - xRender(p) - 0.5) / particleDiameters(p))) * ...
                   (erf( particleDiameters(p) *  (r - yRender(p) + 0.5) / particleDiameters(p)) - erf(particleDiameters(p) * (r - yRender(p) - 0.5) / particleDiameters(p)));
        end
    end    
end

% Flip the image, which contains no noise
IMAGEOUT = flipud(imagePlaceholder);

end
