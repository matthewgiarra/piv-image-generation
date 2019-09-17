function CAMERA_MATRIX = calculate_camera_matrix(Camera)
		
% Camera rotation in world coordinates
rx = Camera.Extrinsic.Rotation.X;
ry = Camera.Extrinsic.Rotation.Y;
rz = Camera.Extrinsic.Rotation.Z;

% Camera position in world coordinates (translation)
tx = Camera.Extrinsic.Translation.X;
ty = Camera.Extrinsic.Translation.Y;
tz = Camera.Extrinsic.Translation.Z;

% Focal length in meters
f = Camera.Intrinsic.FocalLength;

% Number of pixels
image_rows = Camera.Intrinsic.Pixel.Number.Rows;
image_cols = Camera.Intrinsic.Pixel.Number.Columns;

% Pixel sizes in meters
pixel_width_m  = Camera.Intrinsic.Pixel.Width;
pixel_height_m = Camera.Intrinsic.Pixel.Height;

% Rotations about each axis
RX = [1, 0, 0; ...
	  0, cos(rx), -sin(rx); ...
	  0, sin(rx), cos(rx)];

RY = [cos(ry), 0, sin(ry); ...
	  0, 1, 0;  ...
	  -sin(ry), 0, cos(ry)];

RZ = [cos(rz), -sin(rz), 0; ...
	 sin(rz), cos(rz), 0; ...
	 0, 0, 1];

% Total rotation matrix
R = RX * RY * RZ;

% Allocate extrinsic matrix
extrinsic_matrix = zeros(4, 4);

% Populate the rotation part of the matrix
extrinsic_matrix(1:3, 1:3) = R;

% Translation
T = -R * [tx; ty; tz];

% Populate the translation part of the matrix
extrinsic_matrix(1:3, 4) = T;
extrinsic_matrix(4,4) = 1;

% Sensor offsets
xc = image_cols * pixel_width_m / 2;
yc = image_rows * pixel_height_m / 2;

% Pixel size
px = pixel_width_m;
py = pixel_height_m;

% Form the intrinsic matrix
% The sign on the focal length term
% is positive to flip the image into
% image coordinates
intrinsic_matrix = [f / px, 0, xc / px, 0;  ...
					0, -f / py ,yc / py, 0;  ...
					0, 0, 1, 0; ...
					0, 0, 0, 1];

% Calculate the entire camera matrix
% by multiplying the intrinsic and extrinsic
% camera matrices.
CAMERA_MATRIX = intrinsic_matrix * extrinsic_matrix;

end