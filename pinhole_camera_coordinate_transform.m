function [x_cam, y_cam] = pinhole_camera_coordinate_transform(x_world, y_world, z_world, camera_matrix)
% INPUTS
%	x_world = x location of object in world coordinates
%	y_world = '';
%	z_world = '';	
%	camera_matrix = Camera matrix (4 columns x 3 rows)

% OUTPUTS
%	xc = x coordinate of the image in camera coordinate system
%	yc = ''

% Size of position array
array_length = numel(x_world);

% World coordinate matrix
world_coordinates = [(x_world(:))'; (y_world(:))'; (z_world(:))'; ones(1, array_length)];

% Homogeneous camera coordinates
camera_coordinates = camera_matrix * world_coordinates;

% Camera coordinates
x_cam = camera_coordinates(1, :) ./ camera_coordinates(3, :);
y_cam = camera_coordinates(2, :) ./ camera_coordinates(3, :);

end