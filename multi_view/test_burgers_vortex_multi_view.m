function test_burgers_vortex_multi_view(Cameras)

% Seed the number generator
rng(1);

% Number of particles
n_particles = 5E4;

% Particle diameter parameters
particle_diameter_mean = sqrt(8);
particle_diameter_std = sqrt(8) * 0.2;

% Create a normal distribution of particle diameters
particleDiameters = abs(particle_diameter_std * randn(n_particles, 1) ...
    + particle_diameter_mean);

beam_plane_std_dev = 0.02; 
beam_plane_z = 0;

% Time span
tSpan = linspace(0, 1, 11);

% World axis limits in which to generate particles
x_world_limits = 20 * [-1, 1];
y_world_limits = 20 * [-1, 1];
z_world_limits = 0.1 * [-1, 1];

% Count the number of cameras
% Long line because of input checking.
% num_cameras = length(Camera_Parameters);
num_cameras = length(Cameras);

% X world limits (l = lower, u = upper)
xl = x_world_limits(1);
xu = x_world_limits(2);

% Y world limits (l = lower, u = upper)
yl = y_world_limits(1);
yu = y_world_limits(2);

% Z world limits (l = lower, u = upper)
zl = z_world_limits(1);
zu = z_world_limits(2);

% Discrete positions
xo = xl + (xu - xl) * rand(n_particles, 1);
yo = yl + (yu - yl) * rand(n_particles, 1);
zo = zl + (zu - zl) * rand(n_particles, 1);

% Calculate the particle trajectories
[X, Y, Z] = burgersVortex(xo, yo, zo, tSpan);

% Loop over all the time steps
for t = 1 : length(tSpan)

    % [x,y,z] positions at this time point
    x = X(t, :);
    y = Y(t, :);
    z = Z(t, :);
    
    % Calculate particle max intensities from beam profile
    particleMaxIntensities = exp(-(z - beam_plane_z).^2 ./ ...
    (2 * beam_plane_std_dev ^ 2));

    % Inform the user
    fprintf('On frame %d of %d\n', t, length(tSpan));
    
    for k = 1 : num_cameras

        % Grab the current camera matrix
        camera_matrix = Cameras(k).CameraMatrix;

        % Calculate image coordinates
        [x_cam, y_cam] = pinhole_camera_coordinate_transform(x, y, z, camera_matrix);
        n_pixels_rows = Cameras(k).Intrinsic.Pixel.Number.Rows;
        n_pixels_cols = Cameras(k).Intrinsic.Pixel.Number.Columns;
        
        % Sensor gain
        sensorGain = Cameras(k).Intrinsic.Sensor.Gain;
        sensorNoiseStd = Cameras(k).Intrinsic.Sensor.NoiseStd;
        
        % New noise matrix
        noiseMat = abs(sensorNoiseStd * randn(n_pixels_rows, n_pixels_cols));

        particle_image_raw = (...
        generateParticleImage(n_pixels_rows, n_pixels_cols, ...
        x_cam(:), y_cam(:), ...
        particleDiameters, particleMaxIntensities));
  
        % Add the noise
        particle_image_noisy = particle_image_raw + noiseMat;
        
        % Apply sensor gain
        particle_image_uint16 = uint16(sensorGain * double(intmax('uint16')) * particle_image_noisy);
        
        % Make a plot
        subtightplot(2, 2, k, [0.1, 0.1]);
        imagesc(particle_image_uint16);
        axis image;
        set(gca, 'fontsize', 16);
        title(sprintf('Camera %d', k), 'interpreter', 'latex', 'fontsize', 20);
        colormap gray;
        caxis([0, intmax('uint16')]);
        set(gcf, 'color', 'white');
        
        % Make directories if they don't exist
        outDir = Cameras(k).Files.OutDir;
        if(~exist(outDir, 'dir'))
            mkdir(outDir);
        end
        
        % Output path
        outBase = Cameras(k).Files.OutBaseName;
        outPath = fullfile(outDir, sprintf('%s%05d.tiff', outBase, t));
        imwrite(particle_image_uint16, outPath, 'compression', 'none');
        
    end

    % Draw the frame
    drawnow();
end

end




