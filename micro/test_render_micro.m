
% Region size
region_height = 128;
region_width = 128;

% Parse structure: image dimensions in pixels
region_width_pixels  = 128;
region_height_pixels = 128;

% Image magnification (microns per pixel)
pixel_size_microns = 20;

% Optics parameters
% Objective lens
% Magnification (unitless)
objective_magnification = 50;

% Objective lens focal length in microns
focal_length_microns = 4E3;

% Objective lens numerical aperture (unitless)
NA = 0.75;

% Laser wavelength in microns
wavelength_microns = 0.532;

% Fraction above background intensity to render particles
intensity_fraction = 0.0;

% Experiment parameters
% 
% Channel depth in microns
channel_depth_microns = 50;

% Particle diameter in microns
particle_diameter_microns = 1.0;

% Particle concentration (particles per µm^3)
particle_concentration = 5E-3;

% Diffusion
diffusion_std_dev = 3;


% Transform parameters
R = [0, 0, 0];
S = [1, 1, 1];

% Shearing: [sxy, sxz, syx, syz, szx, szy]
sxy = 0;
sxz = 0;
syx = 0;
syz = 0;
szx = 0;
szy = 0;
SH = [sxy, sxz, syx, syz, szx, szy];
T = [0, 0, 0];
Tform = makeAffineTransform_3D(S, R, SH, T);


n_pairs = 100;

% Allocate matrix
imageMatrix1 = zeros(region_height, region_width, n_pairs);
imageMatrix2 = zeros(region_height, region_width, n_pairs);

t1 = tic;

for k = 1 : n_pairs
% Generate image pair.
[imageMatrix1(:, :, k), imageMatrix2(:, :, k)] = ...
    generateImagePair_micro_mc(...
    region_height_pixels, region_width_pixels, particle_diameter_microns ,...
    pixel_size_microns, particle_concentration, channel_depth_microns, ...
    objective_magnification, focal_length_microns, NA, ...
    wavelength_microns, intensity_fraction, diffusion_std_dev, Tform);
end

t2 = toc(t1);

seconds_per_pair = t2 / n_pairs;

fprintf('%0.3f seconds for %d pairs\n', t2, n_pairs);
fprintf('%0.3f seconds per pair\n', seconds_per_pair);


% % Show the images;
% subplot(1, 2, 1);
% imagesc(IMAGE1); axis image; colormap gray;
% c = caxis;
% caxis([0, max(c)]);
% 
% subplot(1, 2, 2);
% imagesc(IMAGE2); axis image; colormap gray;
% caxis([0, max(c)]);


