function CAMERA_PARAMETERS = calculate_camera_parameters(CAMERA_PARAMETERS)

	% Copy the input structure
	Cameras = CAMERA_PARAMETERS.Cameras;

	% Count the number of cameras
	num_cameras = CAMERA_PARAMETERS.NumberOfCameras;
	
	% Target distance
	target_distance = CAMERA_PARAMETERS.TargetPlaneDistance;

	% Loop over all the cameras
	for k = 1 : num_cameras
	
		% Camera rotation in world coordinates
		rx = Cameras(k).Extrinsic.Rotation.X;
		ry = Cameras(k).Extrinsic.Rotation.Y;
		rz = Cameras(k).Extrinsic.Rotation.Z;
	
		% Camera position in world coordinates (translation)
		tx = Cameras(k).Extrinsic.Translation.X;
		ty = Cameras(k).Extrinsic.Translation.Y;
		tz = Cameras(k).Extrinsic.Translation.Z;
	
		% Camera focal length in world units
		focal_length = Cameras(k).Intrinsic.FocalLength;
		
		% Sensor size in world units
		sensor_height_world = Cameras(k).Intrinsic.Sensor.Height;
		sensor_width_world  = Cameras(k).Intrinsic.Sensor.Width;
	
		% Pixel size in world units
		pixel_width_world  = Cameras(k).Intrinsic.Pixel.Width;
		pixel_height_world = Cameras(k).Intrinsic.Pixel.Height;
	
		% Pixel size vector
		pixel_size_vect = [pixel_height_world; pixel_width_world];
	
		% Sensor size vector
		sensor_size_vect = [sensor_height_world; sensor_width_world];
			
		% Calculate camera matrix (pixel coordinates)
		Cameras(k).Camera_Matrix = ...
			calculate_camera_matrix_pixel(rx, ry, rz, tx, ty, tz, focal_length, target_distance, sensor_size_vect, pixel_size_vect);	
		
	end
	
	% Update cameras structure.
	CAMERA_PARAMETERS.Cameras = Cameras;

end