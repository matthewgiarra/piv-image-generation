function test_particle_position_transform_animate(Camera_Arrangement);

% Seed the number generator
rng(1);

% Number of particles
n_particles = 5E3;

% Particle diameter
particle_diameters = 2 * sqrt(8) * ones(n_particles, 1);

% Particle brightness
particle_max_intensity = 1 * ones(n_particles, 1);

% World axis limits
x_world_limits = 0.2 * [-1, 1];
y_world_limits = 0.2 * [-1, 1];
z_world_limits = 0.05 * [-1, 1];

% Count the number of cameras
% Long line because of input checking.
% num_cameras = length(Camera_Parameters);
Camera_Parameters = Camera_Arrangement.Cameras;
num_cameras = length(Camera_Parameters);

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
x = xl + (xu - xl) * rand(n_particles, 1);
y = yl + (yu - yl) * rand(n_particles, 1);
z = zl + (zu - zl) * rand(n_particles, 1);
X = [x'; y'; z'];

%TODO: Y coordinate is flipped when images are generated. 
%   Need to invert Y coordinate in camera matrix.

nframes = 30;
t = linspace(0, 2*pi * (1-1/nframes), nframes);

gap = [0.1, 0.1];
outdir = '~/Desktop/frames_blk_02';

for p = 1 : nframes
    th = 0.5 * sin(t(p));
    rx = [1,0,0; 0, cos(th), -sin(th); 0, sin(th), cos(th)];
    ry = [cos(th), 0, sin(th); 0, 1, 0; -sin(th), 0, cos(th)];
    rz = [cos(th), -sin(th), 0; sin(th), cos(th), 0; 0, 0, 1];
    
    Xp = ry * X;
    
    xp = Xp(1, :);
    yp = Xp(2, :);
    zp = Xp(3, :);
    
    
    for k = 1 : num_cameras

        % Grab the current camera matrix
        camera_matrix = Camera_Parameters(k).Camera_Matrix;

        % Calculate image coordinates
        [x_cam, y_cam] = pinhole_camera_coordinate_transform(xp, yp, zp, camera_matrix);
        n_pixels_rows = Camera_Parameters(k).Intrinsic.Pixel.Number.Rows;
        n_pixels_cols = Camera_Parameters(k).Intrinsic.Pixel.Number.Columns;

        particle_image = (...
        generateParticleImage(n_pixels_rows, n_pixels_cols, ...
        x_cam(:), y_cam(:), ...
        particle_diameters, particle_max_intensity));

        subtightplot(2, 2, k, gap);
        imagesc(particle_image);
        axis image;
    %     set(gca, 'ydir', 'normal');
        set(gca, 'fontsize', 16);
        title(sprintf('Camera %d', k), 'interpreter', 'latex', 'fontsize', 20);

        caxis([0, 3.0]);
        colormap gray

    end
    
    drawnow();
    outname = sprintf('frame_%03d.png', p);
    outpath = fullfile(outdir, outname);
    export_fig('-r150', '-p0.02', outpath);
    
end



% Format axes
% axis image
% set(gca, 'ydir', 'reverse');
% % xlim([0, sensor_width(1) / pixel_size_x(1)]);
% % ylim([0, sensor_height(1) / pixel_size_y(1)]);
% % Release the plot hold
% hold off;



end




