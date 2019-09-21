function [CAMERA_MATRIX, R, t] = getCameraMatrix(Camera)
		
    % Focal length in meters
    f = Camera.FocalLength;

    % Number of pixels
    image_rows = Camera.PixelRows;
    image_cols = Camera.PixelColumns;

    % Pixel sizes in meters
    pixel_width_m  = Camera.PixelWidth;
    pixel_height_m = Camera.PixelHeight;
    
    % Extrinsic matrix
    [extrinsic_matrix, R, t] = cameraLookAtToExtrinsic(Camera.Eye, Camera.Center, Camera.Up);

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
                        0, f / ph ,yc / ph, 0;  ...
                        0, 0, 1, 0; ...
                        0, 0, 0, 1];
                    
    % Calculate the entire camera matrix
    % by multiplying the intrinsic and extrinsic
    % camera matrices.
    CAMERA_MATRIX = intrinsic_matrix * extrinsic_matrix;

end


