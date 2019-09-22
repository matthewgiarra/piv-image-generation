function plotCameraArrangement(varargin)
% Plots the current camera arrangement

% Input parser
p = inputParser;

% Add optional inputs
addParameter(p, 'cameras', defaultCameraArrangement(), @isstruct);
addParameter(p, 'points', [], @isnumeric);

% Parse the arguments
parse(p, varargin{:});

% Results structure
Cameras = p.Results.cameras;
calPoints = p.Results.points;

for n = 1 : length(Cameras)
   
   % Get the camera info
   Camera = Cameras(n);
   [~, R, t] = getCameraMatrix(Camera);
   
   % Convert from camera matrix to orentation matrix
   C = 1 * R\t;
   Rc = inv(R)';
   
   % Plot the camera
   plotCamera('location', C, 'orientation', ...
       Rc, 'label', sprintf('%d', n), 'color', 'w', 'size', 0.05, 'axesvisible', false);
   hold on;
    
end

% Plot the cal target points
if isempty(calPoints)
    [~, ~, ~, x,y,z] = calibrationTarget();
else
    x = calPoints(:, 1);
    y = calPoints(:, 2);
    z = calPoints(:, 3);
end

plot3(x(:),y(:),z(:), '.w', 'markersize', 15, 'markerfacecolor', 'w');
hold off;
axis image;
axis vis3d;
grid on;
box on;

set(gca, 'fontsize', 16);

% Axis labels
xlabel('x (m)', 'fontsize', 16, 'interpreter', 'latex');
ylabel('y (m)', 'fontsize', 16, 'interpreter', 'latex');
zlabel('z (m)', 'fontsize', 16, 'interpreter', 'latex');

% Plot formatting stuff
set(gcf, 'color', 'white');
set(gcf, 'color', 'black');
set(gca, 'color', 'black');
set(gca, 'xcolor', 'white');
set(gca, 'ycolor', 'white');
set(gca, 'zcolor', 'white');


end