function calImages = calTargetMultiView(varargin)

    % Input parser
    p = inputParser;

    % Add optional inputs
    addParameter(p, 'plot', false, @islogical);
    addParameter(p, 'cameras', defaultCameraArrangement(), @isstruct);
    addParameter(p, 'targetOrigin', [0,0,0], @isnumeric);
    
    % Parse the arguments
    parse(p, varargin{:});

    % Results structure
    makePlots = p.Results.plot;

    % Get the cameras
    Cameras = p.Results.cameras;
    
    % Target origin
    targetOrigin = p.Results.targetOrigin;
   
    % Number of targets to render per camera
    nTargs = size(targetOrigin, 1);
    
    for k = 1 : length(Cameras)

        % Ger the current camera
        Camera = Cameras(k);
        
        % Create array for images
        imgArr = zeros(Camera.PixelRows, Camera.PixelColumns, nTargs, 'uint16');
        
        for n = 1 : nTargs

            % Get particle coordinates to render a calibration target
            [x,y,z, xc, yc, zc] = calibrationTarget('origin', targetOrigin(n, :));

            % Create a normal distribution of particle diameters
            particleDiameters = sqrt(8) * ones(numel(x), 1);

            % Calculate particle max intensities from beam profile
            particleMaxIntensities = ones(numel(x), 1);
            
             % Calculate image coordinates
            [x_cam, y_cam] = pinholeTransform(x(:), y(:), z(:), getCameraMatrix(Camera));

            % Render the image and add noise
            particle_image = (...
                generateParticleImage(Camera.PixelRows, Camera.PixelColumns, ...
                  x_cam, y_cam, particleDiameters, particleMaxIntensities)) ...
                  + getSensorNoise(Camera);

            % Apply sensor gain and convert to uint16
            particle_image_uint16 = uint16(Camera.SensorGain * double(intmax('uint16')) * particle_image);

            % Save the image to the array of images
            imgArr(:, :, n) = particle_image_uint16;
            
        end

        % Save results to the output structure
        calImages{k} = imgArr;

        % Plot the images
        if makePlots
        figure(1);
            subtightplot(2, 2, k, [0.1, 0.1]);
            imagesc(imgArr(:, :, 1));
            axis image;
            set(gca, 'fontsize', 16);
            title(sprintf('Camera %d', k), 'interpreter', 'latex', 'fontsize', 20);
            colormap gray;
            caxis([0, intmax('uint16')]);
            set(gca, 'ydir', 'normal');
            set(gcf, 'color', 'white');
        end
    end

    if makePlots
        figure(2);
        plotCameraArrangement(Cameras, 'points', [xc(:), yc(:), zc(:)]);
    end
    
    % Draw the frame
    drawnow();


end




