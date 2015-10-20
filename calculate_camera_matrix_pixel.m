function CAMERA_MATRIX = calculate_camera_matrix_pixel(RX, RY, RZ, TX, TY, TZ, F, SENSOR_SIZE, PIXEL_SIZE);
		
% Allocate extrinsic matrix
extrinsic_matrix = zeros(4, 4);

% Rotations about each axis
RX = [1, 0, 0; ...
	  0, cos(RX), -sin(RX); ...
	  0, sin(RX), cos(RX)];

RY = [cos(RY), 0, sin(RY); ...
	  0, 1, 0;  ...
	  -sin(RY), 0, cos(RY)];

RZ = [cos(RZ), -sin(RZ), 0; ...
	 sin(RZ), cos(RZ), 0; ...
	 0, 0, 1];

% Total rotation matrix
r = RX * RY * RZ;

% Populate the rotation part of the matrix
extrinsic_matrix(1:3, 1:3) = r;

% Populate the translation part of the matrix
extrinsic_matrix(:, 4) = [TX; TY; TZ; 1];

% Sensor offsets
xc = SENSOR_SIZE(2) / 2;
yc = SENSOR_SIZE(1) / 2;

% Pixel size
sx = PIXEL_SIZE(2);
sy = PIXEL_SIZE(1);

% Form the intrinsic matrix
% The sign on the focal length term
% is positive to flip the image into
% image coordinates
intrinsic_matrix = [F / sx, 0, xc / sx, 0;  ...
					0, F / sy ,yc / sy, 0;  ...
					0, 0, 1, 0; ...
					0, 0, 0, 1];

% Calculate the entire camera matrix
% by multiplying the intrinsic and extrinsic
% camera matrices.
CAMERA_MATRIX = intrinsic_matrix * extrinsic_matrix;

end