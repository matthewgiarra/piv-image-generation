function OBJECTIVE_PARAMETERS = load_objective_parameters(objective_name)

% Objective name
OBJECTIVE_PARAMETERS.Name = lower(objective_name);

switch lower(objective_name)
    
    case '10x'       
        % Objective magnification
        OBJECTIVE_PARAMETERS.Magnification = 10;

        % Objective lens focal length in microns
        OBJECTIVE_PARAMETERS.FocalLength = 20E3;

        % Objective lens numerical aperture (unitless)
        OBJECTIVE_PARAMETERS.NA = 0.42;             
    
    case '20x'    
        % Objective magnification
        OBJECTIVE_PARAMETERS.Magnification = 20;

        % Objective lens focal length in microns
        OBJECTIVE_PARAMETERS.FocalLength = 10E3;

        % Objective lens numerical aperture (unitless)
        OBJECTIVE_PARAMETERS.NA = 0.28;        
    case '50x'        
        % Objective magnification
        OBJECTIVE_PARAMETERS.Magnification = 50;

        % Objective lens focal length in microns
        OBJECTIVE_PARAMETERS.FocalLength = 4E3;

        % Objective lens numerical aperture (unitless)
        OBJECTIVE_PARAMETERS.NA = 0.75;
end




end