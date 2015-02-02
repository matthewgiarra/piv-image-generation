%#codegen
function IMAGEOUT = generateParticleImage_3D(HEIGHT, WIDTH, DEPTH, X, Y, Z, PARTICLEDIAMETERS, PARTICLEMAXINTENSITIES)
% This function generate a synthetic 3D PIV particle volume.
% Usage: IMAGEOUT = generateParticleImage_3D(HEIGHT, WIDTH, DEPTH, X, Y, Z, PARTICLEDIAMETERS, PARTICLEMAXINTENSITIES)
%        IMAGEOUT = generateParticleImage_3D_mex(HEIGHT, WIDTH, DEPTH, X, Y, Z, PARTICLEDIAMETERS, PARTICLEMAXINTENSITIES) 
%
% INPUTS:
%   HEIGHT = Number of rows in the volume to generate (integer).
%   WIDTH  = Number of columns in the volume to generate (integer).
%   DEPTH  = Number of depth positions in the volume to generate (integer).
%   X = Column vector containing the non-integer column
%       positions of the particles to generate.
%   Y = Column vector containing the non-integer row
%       positions of the particles to generate.
%   Z = Column vector containing the non-integer depth
%       positions of the particles to generate.
%   PARTICLEDIAMETERS = Column vector containing the non-integer
%       diameters of the particles to generate.
%   PARTICLEMAXINTENSITIES = Column vector containing the maximum
%       intensities of all the particles to generate.
%
% OUTPUTS:
%   IMAGEOUT = [HEIGHT x WIDTH x DEPTH] array containing
%       rendered particle volume.

% Create a placeholder to store the generated image
imagePlaceholder = zeros(HEIGHT, WIDTH, DEPTH); 

% Determine which particles to render. Only render particles whose images lie
% within the image domain. 
particlesToRender = X <= WIDTH & X >= 1 & Y <= HEIGHT & Y >= 1 & Z <= DEPTH & Z >= 1;

% Extract the coordinates of the particles to be rendered.
xRender = X(particlesToRender);
yRender = Y(particlesToRender);
zRender = Z(particlesToRender);

% Extract the max intensities of the particles to be rendered.
particleIntensities = PARTICLEMAXINTENSITIES(particlesToRender);

% Extract the diameters of the particles to be rendered.
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

% Determine the minimum and maximum depth-wise coordinates (foremost and
% backmost voxels) in the image to which each particle
% contributes some intensity.
minRenderedDepth = max(1,   floor(zRender - 0.75 * mean_particle_diameter));
maxRenderedDepth = min(DEPTH, ceil(zRender + 0.75 * mean_particle_diameter));

% Generate the intensities from each particle.
for p = 1 : numberOfParticlesToRender
    for c = minRenderedCols(p) : maxRenderedCols(p)
        for r = minRenderedRows(p) : maxRenderedRows(p)
           for d = minRenderedDepth(p) : maxRenderedDepth(p)
               imagePlaceholder(r, c, d) = imagePlaceholder(r, c) + ...
                particleIntensities(p) * (particleDiameters(p))^2 * pi / 32 *...
                   (erf( particleDiameters(p) *  (c - xRender(p) + 0.5) / particleDiameters(p)) - erf(particleDiameters(p) * (c - xRender(p) - 0.5) / particleDiameters(p))) * ...
                   (erf( particleDiameters(p) *  (r - yRender(p) + 0.5) / particleDiameters(p)) - erf(particleDiameters(p) * (r - yRender(p) - 0.5) / particleDiameters(p))) * ...
                   (erf( particleDiameters(p) *  (d - zRender(p) + 0.5) / particleDiameters(p)) - erf(particleDiameters(p) * (d - zRender(p) - 0.5) / particleDiameters(p)));
           end
        end
    end    
end

% Flip the image, which contains no noise
IMAGEOUT = flipud(imagePlaceholder);

end
