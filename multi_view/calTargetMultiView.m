function calTargetMultiView()

    % Get the cameras
    Cameras = defaultCameraArrangement();

    % Seed the number generator
    rng(1);

    % Number of particles
    particles_per_dot = 5e2;
    
    % Dot spacing
    dL = 0.02;
    
    % Dot diameter in meters
    dotD = 0.005;
 
    % Dot rows and cols
    dotRows = 8;
    dotCols = 8;
    nDots = dotRows * dotCols;
    
    % Width of target (center of first dot to center of last dot)
    targetWidth =  (dotCols - 1) * dL;
    targetHeight = (dotRows - 1) * dL;
    
    % Dot centers
    xv = linspace(-targetWidth/2, targetWidth/2, dotCols);
    yv = linspace(-targetHeight/2, targetHeight/2, dotRows);
    zv = 0;
    
    % Make a grid of dot centers
    [xdots,ydots,zdots] = meshgrid(xv, yv, zv);

    % Create a normal distribution of particle diameters
    particleDiameters = sqrt(8) * ones(particles_per_dot, nDots);
    
    % Calculate particle max intensities from beam profile
    particleMaxIntensities = ones(particles_per_dot, nDots);
    
    % Random angle from dot center
    th = 2 * pi * rand(particles_per_dot, nDots);
    r = dotD/2  * rand(particles_per_dot, nDots);
    zraw = zeros(particles_per_dot, nDots);
    [xraw, yraw, zraw] = pol2cart(th, r, zraw);
    
    x = xraw + (xdots(:))';
    y = yraw + (ydots(:))';
    z = zraw + (zdots(:))';
    
    for k = 1 : length(Cameras)

        % Ger the current camera
        Camera = Cameras(k);
        
        % Camera matrix
        M = getCameraMatrix(Camera);
        
        % Calculate image coordinates
        [x_cam, y_cam] = pinholeTransform(x(:), y(:), z(:), M);
        
        % Render the image and add noise
        particle_image = (...
            generateParticleImage(Camera.PixelRows, Camera.PixelColumns, ...
              x_cam, y_cam, particleDiameters, particleMaxIntensities)) ...
              + getSensorNoise(Camera);

        % Apply sensor gain
        particle_image_uint16 = uint16(Camera.SensorGain * double(intmax('uint16')) * particle_image);

        % Make a plot
        figure(1);
        subtightplot(2, 2, k, [0.1, 0.1]);
        imagesc(particle_image_uint16);
        axis image;
        set(gca, 'fontsize', 16);
        title(sprintf('Camera %d', k), 'interpreter', 'latex', 'fontsize', 20);
        colormap gray;
        caxis([0, intmax('uint16')]);
        set(gca, 'ydir', 'normal');
        set(gcf, 'color', 'white');
        
       
    end

    % TODO: Plot cameras and target in 3D space
    
    % Draw the frame
    drawnow();


end




