function CAMERA_PARAMETERS = default_camera_parameters()

	% Number of cameras
	n_cameras = 4;
    
    % Just some scaling factors
    sr = 0.5;
    st = 0.6;
    
    % Camera positions
    tx = st * [-1, 1, -1, 1];
    ty = st * [1, 1, -1, -1];
    tz = -1 * ones(n_cameras, 1);
    
    % Camera angles
    rx = sr * [-1, -1, 1, 1];
    ry = sr * [-1, 1, -1, 1];
	rz = zeros(n_cameras, 1);
    	
	% Camera focal lengths
	f = 0.028 * ones(n_cameras, 1);

	% Numbers of pixels
	image_rows = 1024 * ones(n_cameras, 1);
	image_cols = 1024 * ones(n_cameras, 1);

	% Pixel sizes (world units)
	pixel_height_world = 1.7E-5 * ones(n_cameras, 1);
	pixel_width_world  = 1.7E-5 * ones(n_cameras, 1);
	
	% Sensor sizes (world units, like mm)
	sensor_width_world  = image_cols .* pixel_width_world;
	sensor_height_world = image_rows .* pixel_height_world;

	% Error check the number of cameras.
	num_cameras = min([n_cameras, length(rx), length(ry), length(rz),...
	 length(tx), length(ty), length(tz), length(f)]);
	 
	 % Save arrangement parameters
	 CAMERA_PARAMETERS.NumberOfCameras = n_cameras; 
 
	 % Loop over cameras
	for k = 1 : num_cameras 
		% Populate  Camera parameters
		CAMERA_PARAMETERS.Cameras(k).Extrinsic.Rotation.X = rx(k);
		CAMERA_PARAMETERS.Cameras(k).Extrinsic.Rotation.Y = ry(k);
		CAMERA_PARAMETERS.Cameras(k).Extrinsic.Rotation.Z = rz(k);

		CAMERA_PARAMETERS.Cameras(k).Extrinsic.Translation.X = tx(k);
		CAMERA_PARAMETERS.Cameras(k).Extrinsic.Translation.Y = ty(k);
		CAMERA_PARAMETERS.Cameras(k).Extrinsic.Translation.Z = tz(k);

		CAMERA_PARAMETERS.Cameras(k).Intrinsic.FocalLength = f(k);
		CAMERA_PARAMETERS.Cameras(k).Intrinsic.Sensor.Width = sensor_width_world(k);
		CAMERA_PARAMETERS.Cameras(k).Intrinsic.Sensor.Height = sensor_height_world(k);
		CAMERA_PARAMETERS.Cameras(k).Intrinsic.Pixel.Height = pixel_height_world(k);
		CAMERA_PARAMETERS.Cameras(k).Intrinsic.Pixel.Width = pixel_width_world(k);
		CAMERA_PARAMETERS.Cameras(k).Intrinsic.Pixel.Number.Rows = image_rows(k);
		CAMERA_PARAMETERS.Cameras(k).Intrinsic.Pixel.Number.Columns = image_cols(k);
	end

	% Update camera parameters with camera matrices.
	CAMERA_PARAMETERS = calculate_camera_parameters(CAMERA_PARAMETERS);

end