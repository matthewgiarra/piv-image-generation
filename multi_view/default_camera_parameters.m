function CAMERAS = default_camera_parameters()

	% Number of cameras
	n_cameras = 4;
    
    % Just some scaling factors
%     st = 0.6;
%     sr = 0.5;
    
    st = 6.0;
    sr = 0.5;
    
    
    % Camera positions
    tx = st * [-1, 1, -1, 1];
    ty = st * [1, 1, -1, -1];
%     tz = -1 * ones(n_cameras, 1);
    tz = -10 * ones(n_cameras, 1);
    
    % Camera angles
    rx = sr * [-1, -1, 1, 1];
    ry = sr * [-1, 1, -1, 1];
	rz = zeros(n_cameras, 1);
    
    % gain = 1 means that a floating point
    % pixel value of 1 corresponds to a saturated pixel.
    % gain = 0 will make the whole image black.
    % You have to figure out a good value here by trial and error.
    defaultSensorGain = 0.6;
    defaultSensorNoiseStd = 0.05;
    	
	% Camera focal lengths
	f = 0.028 * ones(n_cameras, 1);

	% Numbers of pixels
	image_rows = 1024 * ones(n_cameras, 1);
	image_cols = 1024 * ones(n_cameras, 1);

	% Pixel sizes (world units)
    % These values are for a Photron camera, but I forget which one.
	pixel_height_world = 1.7E-5 * ones(n_cameras, 1);
	pixel_width_world  = 1.7E-5 * ones(n_cameras, 1);
	
    % Sensor gain
    sensorGain = defaultSensorGain * ones(n_cameras, 1);
    sensorNoiseStd = defaultSensorNoiseStd * ones(n_cameras, 1);
    
   	% Error check the number of cameras.
	num_cameras = min([n_cameras, length(rx), length(ry), length(rz),...
	 length(tx), length(ty), length(tz), length(f)]);
	      
     outBaseDir = 'images';
     outBaseName = 'frame_';
      
	 % Loop over cameras
	for k = 1 : num_cameras 
		% Populate  Camera parameters
		CAMERAS(k).Extrinsic.Rotation.X = rx(k);
		CAMERAS(k).Extrinsic.Rotation.Y = ry(k);
		CAMERAS(k).Extrinsic.Rotation.Z = rz(k);

		CAMERAS(k).Extrinsic.Translation.X = tx(k);
		CAMERAS(k).Extrinsic.Translation.Y = ty(k);
		CAMERAS(k).Extrinsic.Translation.Z = tz(k);

		CAMERAS(k).Intrinsic.FocalLength = f(k);
		CAMERAS(k).Intrinsic.Pixel.Height = pixel_height_world(k);
		CAMERAS(k).Intrinsic.Pixel.Width = pixel_width_world(k);
		CAMERAS(k).Intrinsic.Pixel.Number.Rows = image_rows(k);
		CAMERAS(k).Intrinsic.Pixel.Number.Columns = image_cols(k);
        CAMERAS(k).Intrinsic.Sensor.NoiseStd = sensorNoiseStd(k);
        CAMERAS(k).Intrinsic.Sensor.Gain = sensorGain(k);
        CAMERAS(k).CameraMatrix = calculate_camera_matrix(CAMERAS(k));
        
        % Camera save locations
        CAMERAS(k).Files.OutDir = fullfile(outBaseDir, sprintf('cam%02d', k));
        CAMERAS(k).Files.OutBaseName = outBaseName;
        
	end


end