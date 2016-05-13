function diffusion_std_dev_pix = calculate_particle_diffusion(...
    channel_width_microns, channel_depth_microns, ...
    flow_rate_ul_min, objective_magnification, ...
    pixel_size_microns, T_kelvin, dp_microns, viscosity_pas, dx_target_pix);

% % Channel dimensions
% channel_width_microns = 5E3;
% channel_depth_microns = 100;

% % Flow rate in µL/min
% flow_rate_ul_min = 5;

% % Objective magnification 
% objective_magnification = 60;

% % Pixel size in microns
% pixel_size_microns = 10;

% % Temperature in kelvin
% T_kelvin = 300;

% % Particle diameter in microns
% dp_microns = 0.1;

% % Viscosity in Pa*S
% viscosity_pas = 1.12E-3;

% % Desired pixel displacement in pixels per frame
% dx_target_pix = 15;

% Calculate channel cross section in microns^2
cross_section_microns = channel_width_microns * channel_depth_microns;

% Convert channel cross section to meters^2
cross_section_meters = cross_section_microns * 10^-12;

% Convert pixel size to meters
pixel_size_meters = pixel_size_microns * 10^-6;

% Convert flow rate to m^3 / sec
flow_rate_meters_sec = flow_rate_ul_min / 60 * 10^-9;

% Observed flow velocity in meters per second
velocity_meters_sec = flow_rate_meters_sec / cross_section_meters ...
    * objective_magnification;

% Average flow velocity in image coordinates
velocity_pixels_sec = 1 * velocity_meters_sec / pixel_size_meters;

% DT to achieve target displacement
dt_sec = dx_target_pix / velocity_pixels_sec;

% Diffusion constant in m^2 / sec
diffusion_constant_meters_sec = calculate_diffusion_constant(...
    T_kelvin, dp_microns, viscosity_pas);

% Diffusion standard deviation
diffusion_std_dev_pix = sqrt(2 * diffusion_constant_meters_sec * dt_sec) ...
    * objective_magnification / pixel_size_meters;

end



