function makeImages(varargin)

% Seed the number generator
rng(1);

% Input parser
p = inputParser;

% Add optional inputs
addParameter(p, 'cameras', defaultCameraArrangement(), @isstruct);
addParameter(p, 'outdir', '.', @isstr);
addParameter(p, 'outbase', 'frame_', @isstr);
addParameter(p, 'extension', 'tiff', @isstr);
addParameter(p, 'zeros', 4, @isnumeric);
addParameter(p, 'velocityFunction', @burgersVortex, @isFunctionHandle);
addParameter(p, 'xrange', [-1, 1]);
addParameter(p, 'yrange', [-1, 1]);
addParameter(p, 'zrange', [-0.1, 0.1]);
addParameter(p, 'particleConcentration', 2.5e5, @isnumeric);
addParameter(p, 'tspan', linspace(0,0.01, 20), @isnumeric);
addParameter(p, 'particleDiameterMean', 1.5*sqrt(8), @isnumeric);
addParameter(p, 'particleDiameterStdDev', 0.15 * sqrt(8), @isnumeric);
addParameter(p, 'beamStdDev', 0.05, @isnumeric);
addParameter(p, 'BeamPlaneZ', 0, @isnumeric);
addParameter(p, 'plot', false, @islogical);

% Parse the arguments
parse(p, varargin{:});

Cameras = p.Results.cameras;
out_root = p.Results.outdir;
out_base = p.Results.outbase;
out_ext = p.Results.extension;
nZeros = p.Results.zeros;
velocityFunction = p.Results.velocityFunction;
xrange = p.Results.xrange;
yrange = p.Results.yrange;
zrange = p.Results.zrange;
particle_concentration = p.Results.particleConcentration;
tSpan = p.Results.tspan;
particle_diameter_mean = p.Results.particleDiameterMean;
particle_diameter_std = p.Results.particleDiameterStdDev;
beam_plane_std_dev = p.Results.beamStdDev;
makePlots = p.Results.plot;
beam_plane_z = p.Results.BeamPlaneZ;

% Format string
fmtStr = sprintf('%%0%dd', nZeros);
outNameFmt = sprintf('%s%s.%s', out_base, fmtStr, out_ext);

testVolume = diff(xrange) * diff(yrange) * diff(zrange);
n_particles = particle_concentration * testVolume;

% Create a normal distribution of particle diameters
particleDiameters = abs(particle_diameter_std * randn(n_particles, 1) ...
    + particle_diameter_mean);

% Discrete positions
xo = xrange(1) + (xrange(2) - xrange(1)) * rand(n_particles, 1);
yo = yrange(1) + (yrange(2) - yrange(1)) * rand(n_particles, 1);
zo = zrange(1) + (zrange(2) - zrange(1)) * rand(n_particles, 1);

% Calculate the particle trajectories
[X, Y, Z] = velocityFunction(xo, yo, zo, tSpan);

% Count the number of cameras
num_cameras = length(Cameras);

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
    
    % Loop over cameras
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
        
        if makePlots        
            % Make a plot
            subtightplot(2, 2, k, [0.1, 0.1]);
            imagesc(particle_image_uint16);
            axis image;
            set(gca, 'fontsize', 16);
            title(sprintf('Camera %d', k), 'interpreter', 'latex', 'fontsize', 20);
            colormap gray;
            caxis([0, intmax('uint16')]);
            set(gcf, 'color', 'white');
        end
        
        % Output path
        out_dir = fullfile(out_root, sprintf('Cam%d', k));
        if(~exist(out_dir, 'dir'))
            mkdir(out_dir);
        end
        
        out_path = fullfile(out_dir, sprintf(outNameFmt, t));
        imwrite(particle_image_uint16, out_path);
        
    end
    
    % Draw the frame
    drawnow();
end

end




