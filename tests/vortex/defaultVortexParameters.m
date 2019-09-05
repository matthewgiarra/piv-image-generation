function VORTEXPARAMETERS = defaultVortexParameters(VORTEXMODEL);

isRankine = ~isempty(regexpi(VORTEXMODEL, 'ran'));
isLamb = ~isempty(regexpi(VORTEXMODEL, 'lam'));

if isRankine
% Rankine vortex parameters
    VORTEXPARAMETERS.Vorticity = 1.00;
    VORTEXPARAMETERS.CoreRadius = 50;
    VORTEXPARAMETERS.VortexRadius = 100;
    VORTEXPARAMETERS.XC = 0.500;
    VORTEXPARAMETERS.YC = 0.500;
    VORTEXPARAMETERS.Angle = 0;
    VORTEXPARAMETERS.Circulation = 250;
    VORTEXPARAMETERS.Viscosity = 10;

elseif isLamb
% Lamb-Oseen vortex parameters
    VORTEXPARAMETERS.PeakVelocity = 1.00;
    VORTEXPARAMETERS.CoreRadius = 25;
    VORTEXPARAMETERS.VortexRadius = 200;
    VORTEXPARAMETERS.XC = 500;
    VORTEXPARAMETERS.YC = 500;
    VORTEXPARAMETERS.Angle = 0;
    VORTEXPARAMETERS.PropagationVelocity = 1;
    
else
    VORTEXPARAMETERS = '';
    disp('Invalid vortex model type specified.');
    disp('Specify either "rankine" or "lamb". ');
    
end



end