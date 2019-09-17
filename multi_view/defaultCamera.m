function CAMERA = defaultCamera()

    % Instantiate an empty structure to
    % hold the camera parameters
    CAMERA = struct();
    
    % Translation along [x,y,z] axes in meters
    % The vector [tx, ty, tz] specifies the vector
    % pointing from the origin of the world coordinate system
    % to the origin of the camera coordinate system.
    tx = 0; ty = 0; tz = 0;
    
    % Rotation about [x,y,z] axes in radians
    % The vector [rx, ry, rz] specifies rotations
    % about the axes of the camera coordinate system,
    % not the world coordinate system.
    rx = 0; ry = 0; rz = 0;
    
    % Populate  Camera parameters
    CAMERA.FocalLength      = 0.028;  % Meters
    CAMERA.PixelHeight      = 1.7E-5; % Meters
    CAMERA.PixelWidth       = 1.7E-5; % Meters
    CAMERA.PixelRows        = 1024;   % Pixels
    CAMERA.PixelColumns     = 1024;   % Pixels
    CAMERA.SensorNoiseStd   = 0.05;   % Unitless
    CAMERA.SensorGain       = 0.6;    % Unitless
    CAMERA.Translation      = [tx, ty, tz]; % Meters
    CAMERA.Rotation         = [rx, ry, rz]; % Radians
      
end



