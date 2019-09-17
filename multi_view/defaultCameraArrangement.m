function CAMERAS = defaultCameraArrangement()
    
    % Camera positions
    tx = 6.0 * [-1, 1, -1, 1];
    ty = 6.0 * [1,  1, -1, -1];
    tz = [-10, -10, -10, -10];
    
    % Camera angles
    rx = 0.5 * [-1, -1, 1, 1];
    ry = 0.5 * [-1, 1, -1, 1];
	rz = [0, 0, 0, 0];
    
    % Number of cameras
    nCameras = length(tx);
    
	 % Loop over cameras
    for k = 1 : nCameras 
        
        % Instantiate a camera 
        % with default parameters
        C = defaultCamera;
        
		% Populate Camera parameters
        C.Translation = [tx(k), ty(k), tz(k)];
		C.Rotation    = [rx(k), ry(k), rz(k)];
	  
        % Append to the structure
        CAMERAS(k) = C;
        
    end
end