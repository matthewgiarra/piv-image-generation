function testBurgersVortexMultiView()

% Seed the number generator
rng(1);

% Cameras
Cameras = defaultCameraArrangement();

% Output directory 
out_root = 'images';

% Number of particles
% n_particles = 10E4;
particle_concentration = 2.5E5;

% Particle diameter parameters
particle_diameter_mean = 1.5*sqrt(8);
particle_diameter_std  = 0.1 * particle_diameter_mean;

% World axis limits in which to generate particles
x_world_limits = 1 * [-1, 1];
y_world_limits = 1 * [-1, 1];
z_world_limits = 0.1 * [-1, 1];

testVolume = diff(x_world_limits) * diff(y_world_limits) * diff(z_world_limits);
n_particles = particle_concentration * testVolume;

% Create a normal distribution of particle diameters
particleDiameters = abs(particle_diameter_std * randn(n_particles, 1) ...
    + particle_diameter_mean);

beam_plane_std_dev = 0.05; 
beam_plane_z = 0;

% Time span
tSpan = linspace(0, 0.05, 2);

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
    % in world coordinates (lab frame, meters)
    x = X(t, :);
    y = Y(t, :);
    z = Z(t, :);
    
    % Calculate particle max intensities from beam profile
    particleMaxIntensities = exp(-(z - beam_plane_z).^2 ./ ...
    (2 * beam_plane_std_dev ^ 2));

    % Inform the user
    fprintf('On frame %d of %d\n', t, length(tSpan));
    
    for k = 1 : num_cameras

        % Ger the current camera
        Camera = Cameras(k);
        
        % Calculate image coordinates
        [x_cam, y_cam] = pinholeTransform(x, y, z, getCameraMatrix(Camera));
        
        % Render the image and add noise
        particle_image = (...
            generateParticleImage(Camera.PixelRows, Camera.PixelColumns, ...
              x_cam, y_cam, particleDiameters, particleMaxIntensities)) ...
              + getSensorNoise(Camera);
        
        % Apply sensor gain
        particle_image_uint16 = uint16(Camera.SensorGain * double(intmax('uint16')) * particle_image);
        
        % Make a plot
        subtightplot(2, 2, k, [0.1, 0.1]);
        imagesc(particle_image_uint16);
        axis image;
        set(gca, 'fontsize', 16);
        title(sprintf('Camera %d', k), 'interpreter', 'latex', 'fontsize', 20);
        colormap gray;
        caxis([0, intmax('uint16')]);
        set(gcf, 'color', 'white');
        
        % Output path
        out_dir = fullfile(out_root, sprintf('Cam%d', k));
        if(~exist(out_dir, 'dir'))
            mkdir(out_dir);
        end
        out_path = fullfile(out_dir, sprintf('cam%d_frame_%05d.tiff', k, t));
%         imwrite(particle_image_uint16, out_path, 'compression', 'none');
        
    end

    % Draw the frame
    drawnow();
end

end




