function CAMERAS = defaultCameraArrangement()
   
    % Camera positions
    tx = 0.5 * [-1, 1, -1, 1];
    ty = 0.5 * [1,  1, -1, -1];
    tz = 1 * [1, 1, 1, 1];
    
    % Number of cameras
    nCameras = length(tx);
    
	 % Loop over cameras
    for k = 1 : nCameras 
        
        % Instantiate a camera 
        % with default parameters
        C = defaultCamera;
        
		% Populate Camera parameters
        C.Eye = [tx(k), ty(k), tz(k)];
		
        % Append to the structure
        CAMERAS(k) = C;
    end
end