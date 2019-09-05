%#codegen
function IMAGE_OUT = generateParticleImage(HEIGHT, WIDTH, X, Y, ...
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

% Square root of 8; just calculate this once
sqrt8 = sqrt(8);

% Define render fraction
% This is a multiple of the particle
% diameter that specifies how far away from the 
% particle center to render.
render_fraction = 0.75;

% Create a placeholder to store the generated image
IMAGE_OUT = zeros(HEIGHT, WIDTH); 

% Only render particles whose intensities are above some threshold,
% which we specify here as the value of a zero-mean Gaussian function
% at two standard deviations (so that the "width" of a symmetric 
% Gaussian beam is eight times the beam's standard deviation).
% cutoff_intensity = exp(-6);

% Determine the miniumum and maximum columns (leftmost and rightmost pixels) in the image
% to which each particle contributes some intensity,
% fractional values
minRenderedCols = floor(X - render_fraction * PARTICLE_DIAMETERS);
maxRenderedCols =  ceil(X + render_fraction * PARTICLE_DIAMETERS);

% Determine the minimum and maximum rows (topmost and bottommost pixels) in
% the image to which each particle contributes some intensity, 
% fractional values
minRenderedRows = floor(Y - render_fraction * PARTICLE_DIAMETERS);
maxRenderedRows =  ceil(Y + render_fraction * PARTICLE_DIAMETERS);

% Determine which particles contribute intensity to the image
render_particle = minRenderedCols <= WIDTH & ...
               	  maxRenderedCols >= 1 & ...
                  minRenderedRows <= HEIGHT & ...
                  maxRenderedRows >= 1;

% Incices to render
render_inds = find(render_particle);
              
% Number of particles to render
num_to_render = length(render_inds);

% Generate the intensities from each particle.
for p = 1 : num_to_render
    
    % Particle index
    ind = render_inds(p);
 
    % Inform the user
    fprintf(1, 'On particle %d of %d\n', p, num_to_render);
    
    % Loop over all the pixels to which the particle contributes intensity
    for c = minRenderedCols(ind) : maxRenderedCols(ind)
        for r = minRenderedRows(ind) : maxRenderedRows(ind)
            
            % Radius from the particle center
            render_radius = sqrt( (c - X(ind))^2 + (r - Y(ind)^2));
            
            % Boolean for whether to render the particle
            render_pixel = c >= 1 ...
                && c <= WIDTH ...
                && r >= 1 ...
                && r <= HEIGHT...        
                && render_radius < render_fraction * PARTICLE_DIAMETERS(ind);
           
            % Render the pixel if it meets the criteria
            if render_pixel              
                IMAGE_OUT(r, c) = IMAGE_OUT(r, c) + ...
                   PARTICLE_MAX_INTENSITIES(ind) * (PARTICLE_DIAMETERS(ind))^2 * pi / 32 *...
                       ( erf( sqrt8 *  (c - X(ind) - 0.5) ...
                       / PARTICLE_DIAMETERS(ind) ) - erf(sqrt8 * ...
                       (c - X(ind) + 0.5) / PARTICLE_DIAMETERS(ind))) * ...
                       ...
                       (erf( sqrt8 *  (r - Y(ind) - 0.5) ...
                       / PARTICLE_DIAMETERS(ind)) - erf(sqrt8 * ...
                       (r - Y(ind) + 0.5) / PARTICLE_DIAMETERS(ind)));               
            end      
        end
    end    
end

% % Keep the image from over-saturating by limiting its maximum value to one.
% if max(IMAGE_OUT(:)) > 1
%     IMAGE_OUT = IMAGE_OUT ./ max(IMAGE_OUT(:));
% end

% % Flip the image, which contains no noise
% IMAGEOUT = flipud(imagePlaceholder);

end
