function CAMERA_PARAMETERS = default_camera_parameters();

	% Number of cameras
	n_cameras = 3;

	% Distance from camera arrangement centroid to the target plane.
	L = 1;

	% Camera angles
	rx = 0 * ones(n_cameras, 1);
	ry = 0 * ones(n_cameras, 1);
	rz = 0 * ones(n_cameras, 1);

	% Camera positions
	tx = 0.1 * [0, -1, 1];
	ty = 0.1 * [1, -1, -1];
	tz = L   * [1, 1, 1];

	% Camera focal lengths
	f = 0.028 * ones(n_cameras, 1);

	% Sensor sizes (world units, like mm)
	sensor_width_world  = 1.74 * 1E-2 * ones(n_cameras, 1);;
	sensor_height_world = 1.74 * 1E-2 * ones(n_cameras, 1);;

	% Pixel sizes (world units)
	pixel_height_world = 1.7E-5 * ones(n_cameras, 1);;
	pixel_width_world = 1.7E-5 * ones(n_cameras, 1);;

	% Numbers of pixels
	image_rows = 512 * ones(n_cameras, 1);;
	image_cols = 512 * ones(n_cameras, 1);;

	% Error check the number of cameras.
	num_cameras = min([n_cameras, length(rx), length(ry), length(rz),...
	 length(tx), length(ty), length(tz), length(f)]);
	 
	 % Save arrangement parameters
	 CAMERA_PARAMETERS.NumberOfCameras = n_cameras;
	 CAMERA_PARAMETERS.TargetPlaneDistance = L; 
 
	 % Loop over cameras
	for k = 1 : num_cameras 
		% Populate  Camera parameters
		CAMERA_PARAMETERS.Cameras(k).Extrinsic.Rotation.X = rx(k);
		CAMERA_PARAMETERS.Cameras(k).Extrinsic.Rotation.Y = ry(k);
		CAMERA_PARAMETERS.Cameras(k).Extrinsic.Rotation.Z = ry(k);

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