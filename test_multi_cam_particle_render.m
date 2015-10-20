function test_multi_cam_particle_render();

% Number of cameras
n_cameras = 3;

% Camera angles
rx = 0 * [1, 1, 1];
ry = 0 * [1, 1, 1];
rz = 0 * [1, 1, 1];

% Camera positions
tx = 0.005 * [0, -1, 1]; 
ty = 0.010 * [0, -1, -1];
tz = 2.0 * [1, 1, 1];

% Camera focal lengths
f = 0.028 * [1, 1, 1];

% Sensor sizes (mm)
sensor_width  = 1.74 * 1E-2 * [1, 1, 1];
sensor_height = 1.74 * 1E-2 * [1, 1, 1];

% Pixel sizes
pixel_height = 1.7E-5 * [1, 1, 1];
pixel_width = 1.7E-5 * [1, 1, 1];

% Error check the number of cameras.
num_cameras = min([n_cameras, length(rx), length(ry), length(rz),...
 length(tx), length(ty), length(tz), length(f)]);
 
for k = 1 : num_cameras 

	% Camera parameters
	Parameters.Cameras(k).Extrinsic.Rotation.X = rx(k);
	Parameters.Cameras(k).Extrinsic.Rotation.Y = ry(k);
	Parameters.Cameras(k).Extrinsic.Rotation.Z = ry(k);

	Parameters.Cameras(k).Extrinsic.Translation.X = tx(k);
	Parameters.Cameras(k).Extrinsic.Translation.Y = ty(k);
	Parameters.Cameras(k).Extrinsic.Translation.Z = tz(k);

	Parameters.Cameras(k).Intrinsic.FocalLength = f(k);
	Parameters.Cameras(k).Intrinsic.Sensor.Width = sensor_width(k);
	Parameters.Cameras(k).Intrinsic.Sensor.Height = sensor_width(k);
	Parameters.Cameras(k).Intrinsic.Pixel.Height = pixel_height(k);
	Parameters.Cameras(k).Intrinsic.Pixel.Width = pixel_width(k);
end

% Calculate and display particle positions
test_particle_position_transform(Parameters.Cameras);


end