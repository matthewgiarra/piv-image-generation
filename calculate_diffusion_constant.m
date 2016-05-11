function diffusion_constant = calculate_diffusion_constant(T_kelvin,...
    dp_microns, viscosity_pas)
% Calculate the diffusion coefficient in microns^2 / sec

% Boltzman's constant
% m^2 * kg / (s^2 * K)
kb = 1.38064852E-23; 

% Particle diameter in meters
dp_meters = dp_microns / 10^6;

% Calculate diffusion (m^2 / s)
% The factor of 10^12 converts meters^2 to microns^2
% In the original code this was 10^6 not 10^12, which
% seems wrong to me. 
diffusion_constant = kb * T_kelvin ./ ...
    (6 * pi * viscosity_pas .* dp_meters / 2);

end