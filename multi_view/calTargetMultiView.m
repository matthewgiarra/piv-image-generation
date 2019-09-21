function calTargetMultiView()

    % Get the cameras
    Cameras = defaultCameraArrangement();
    
    % Get particle coordinates to render a calibration target
    [x,y,z] = calibrationTarget();
  
    % Create a normal distribution of particle diameters
    particleDiameters = sqrt(8) * ones(numel(x), 1);
    
    % Calculate particle max intensities from beam profile
    particleMaxIntensities = ones(numel(x), 1);
  
    for k = 1 : length(Cameras)

        % Ger the current camera
        Camera = Cameras(k);
                
        % Calculate image coordinates
        [x_cam, y_cam] = pinholeTransform(x(:), y(:), z(:), getCameraMatrix(Camera));
        
        % Render the image and add noise
        particle_image = (...
            generateParticleImage(Camera.PixelRows, Camera.PixelColumns, ...
              x_cam, y_cam, particleDiameters, particleMaxIntensities)) ...
              + getSensorNoise(Camera);

        % Apply sensor gain and convert to uint16
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




