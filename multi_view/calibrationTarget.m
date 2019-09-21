function [x,y,z] = calibrationTarget(varargin)
% Makes [x,y,z] coordinates that can be passed to generateParticleImage 
% to create an image of a camera calibration target.

    % Input parser
    p = inputParser;

    % Add optional inputs
    addParameter(p, 'dotSpacing', 0.0254, @isnumeric);
    addParameter(p, 'dotDiameter', 0.005, @isnumeric);
    addParameter(p, 'rows', 9, @isnumeric);
    addParameter(p, 'columns', 9, @isnumeric);
    addParameter(p, 'origin', [0,0,0], @isnumeric);
    addParameter(p, 'particlesPerDot', 1e3, @isnumeric);
    
    % Parse the arguments
    parse(p, varargin{:});

    % Results structure
    dot_spacing_m = p.Results.dotSpacing;
    dot_diameter_m   = p.Results.dotDiameter;
    dot_rows = p.Results.rows;
    dot_cols = p.Results.columns;
    target_origin   = p.Results.origin;
    particles_per_dot = p.Results.particlesPerDot;
 
    % Number of dots
    nDots = dot_rows * dot_cols;
    
    % Width of target (center of first dot to center of last dot)
    targetWidth =  (dot_cols - 1) * dot_spacing_m;
    targetHeight = (dot_rows - 1) * dot_spacing_m;
    
    % Dot centers
    xv = linspace(-targetWidth/2, targetWidth/2, dot_cols)   + target_origin(1);
    yv = linspace(-targetHeight/2, targetHeight/2, dot_rows) + target_origin(2);
    zv = target_origin(3);
    
    % Make a grid of dot centers
    [xdots, ydots, zdots] = meshgrid(xv, yv, zv);
    
    % Reshape the dot center arrays into vectors
    xdotsVect = xdots(:);
    ydotsVect = ydots(:);
    zdotsVect = zdots(:);
    
    % Allocate array to hold all the [x,y,z] points
    x = zeros(particles_per_dot, nDots);
    y = zeros(particles_per_dot, nDots);
    z = zeros(particles_per_dot, nDots);
    
    % Make all the dots
    for n = 1 : nDots
        
        % Random angle from dot center
        th = 2 * pi * rand(particles_per_dot, 1);
        r = dot_diameter_m/2  * rand(particles_per_dot, 1);
        zraw = zeros(particles_per_dot, 1);
        [xraw, yraw, ~] = pol2cart(th, r, zraw);
        x(:, n) = xraw + xdotsVect(n);
        y(:, n) = yraw + ydotsVect(n);
        z(:, n) = zraw + zdotsVect(n);
    end
end
