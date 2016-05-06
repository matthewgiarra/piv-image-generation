function OBJECTIVE_PARAMETERS = load_objective_parameters(objective_name)

% Objective name
OBJECTIVE_PARAMETERS.Name = lower(objective_name);

switch lower(objective_name)
    
    case '10x'
        
        % 10X Mitutoyo Plan Apo Infinity Corrected Long WD Objective

        % Objective magnification
        OBJECTIVE_PARAMETERS.Magnification = 10;
        
        % Objective lens numerical aperture (unitless)
        OBJECTIVE_PARAMETERS.NA = 0.28;
        
        % Objective lens working distance (microns)
        OBJECTIVE_PARAMETERS.WorkingDistance = 33.5E3;

        % Objective lens focal length in microns
        OBJECTIVE_PARAMETERS.FocalLength = 20E3;
 
    case '20x'
        
        % 20X Mitutoyo Plan Apo Infinity Corrected Long WD Objective         

        % Objective magnification
        OBJECTIVE_PARAMETERS.Magnification = 20;
        
        % Objective lens numerical aperture (unitless)
        OBJECTIVE_PARAMETERS.NA = 0.42;  
        
        % Objective lens working distance (microns)
        OBJECTIVE_PARAMETERS.WorkingDistance = 20.0E3;

        % Objective lens focal length in microns
        OBJECTIVE_PARAMETERS.FocalLength = 10E3;
              
    case '50x'
        
        % Objective magnification
        OBJECTIVE_PARAMETERS.Magnification = 50;
        
        % Objective lens numerical aperture (unitless)
        OBJECTIVE_PARAMETERS.NA = 0.55;
        
        % Objective lens working distance (microns)
        OBJECTIVE_PARAMETERS.WorkingDistance = 13.0E3;

        % Objective lens focal length in microns
        OBJECTIVE_PARAMETERS.FocalLength = 4E3;

end




end