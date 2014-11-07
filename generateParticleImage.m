%#codegen
function IMAGEOUT = generateParticleImage(HEIGHT, WIDTH, X, Y, PARTICLEDIAMETER, PARTICLEMAXINTENSITIES)

% Create a placeholder to store the generated image
imagePlaceholder = zeros(HEIGHT, WIDTH); 

% Determine which particles to render. Only render particles whose images lie
% within the image domain. 
particlesToRender = X <= WIDTH & X >= 1 & Y <= HEIGHT & Y >= 1;

% Extract the coordinates and brightness of the particles to be rendered.
xRender = X(particlesToRender);
yRender = Y(particlesToRender);
particleIntensities = PARTICLEMAXINTENSITIES(particlesToRender);

% Determine the miniumum and maximum columns (leftmost and rightmost pixels) in the image
% to which each particle contributes some intensity.
minRenderedCols = max(1, floor(xRender - 0.75 * PARTICLEDIAMETER));
maxRenderedCols = min(WIDTH, ceil(xRender + 0.75 * PARTICLEDIAMETER));

% Determine the minimum and maximum rows (topmost and bottommost pixels) in
% the image to which each particle contributes some intensity.
minRenderedRows = max(1, floor(yRender - 0.75 * PARTICLEDIAMETER));
maxRenderedRows = min(HEIGHT, ceil(yRender + 0.75 * PARTICLEDIAMETER));

% Determine the number of particles to render
numberOfParticlesToRender = length(xRender);

% Generate the intensities from each particle.
for p = 1 : numberOfParticlesToRender
    for c = minRenderedCols(p) : maxRenderedCols(p)
        for r = minRenderedRows(p) : maxRenderedRows(p)
            imagePlaceholder(r, c) = imagePlaceholder(r, c) + ...
                particleIntensities(p) * PARTICLEDIAMETER^2 * pi / 32 *...
                   (erf( PARTICLEDIAMETER *  (c - xRender(p) + 0.5) / PARTICLEDIAMETER) - erf(PARTICLEDIAMETER * (c - xRender(p) - 0.5) / PARTICLEDIAMETER)) * ...
                   (erf( PARTICLEDIAMETER *  (r - yRender(p) + 0.5) / PARTICLEDIAMETER) - erf(PARTICLEDIAMETER * (r - yRender(p) - 0.5) / PARTICLEDIAMETER));
        end
    end    
end

% Flip the image, which contains no noise
IMAGEOUT = flipud(imagePlaceholder);

end
