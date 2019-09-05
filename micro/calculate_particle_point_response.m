function point_response = calculate_particle_point_response(...
    MAGNIFICATION, NA, WAVELENGTH)
% Equation from Olsen and Adrian 2000


% f# number
f = 2 * NA;

% Point response of diffraction limited optics
% microscope model. Equation 2, Olsen & Adrian 2000
point_response = 2.44 * (MAGNIFICATION + 1) * WAVELENGTH * f;


end