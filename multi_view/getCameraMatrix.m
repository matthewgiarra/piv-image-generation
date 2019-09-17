function CAMERA_MATRIX = getCameraMatrix(Camera)
		
    % Camera rotation in world coordinates
    rx = Camera.Rotation(1);
    ry = Camera.Rotation(2);
    rz = Camera.Rotation(3);

    % Camera position in world coordinates (translation)
    tx = Camera.Translation(1);
    ty = Camera.Translation(2);
    tz = Camera.Translation(3);

    % Focal length in meters
    f = Camera.FocalLength;

    % Number of pixels
    image_rows = Camera.PixelRows;
    image_cols = Camera.PixelColumns;

    % Pixel sizes in meters
    pixel_width_m  = Camera.PixelWidth;
    pixel_height_m = Camera.PixelHeight;

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

    % Translation matrix
    T = -R * [tx; ty; tz];

    % Populate the translation part of the matrix
    extrinsic_matrix(1:3, 4) = T;
    extrinsic_matrix(4,4) = 1;

    % Sensor offsets
    xc = image_cols * pixel_width_m / 2;
    yc = image_rows * pixel_height_m / 2;

    % Pixel size
    % Just renaming so the code is legible
    pw = pixel_width_m;
    ph = pixel_height_m;

    % Form the intrinsic matrix
    % The sign on the focal length term
    % is positive to flip the image into
    % image coordinates
    intrinsic_matrix = [f / pw, 0, xc / pw, 0;  ...
                        0, -f / ph ,yc / ph, 0;  ...
                        0, 0, 1, 0; ...
                        0, 0, 0, 1];

    % Calculate the entire camera matrix
    % by multiplying the intrinsic and extrinsic
    % camera matrices.
    CAMERA_MATRIX = intrinsic_matrix * extrinsic_matrix;

end