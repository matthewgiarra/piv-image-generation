function [x_world, y_world] = calculate_camera_field_of_view(CameraParameters, Z, L);

% Number of cameras
tx = CameraParameters.Extrinsic.Translation.X;
ty = CameraParameters.Extrinsic.Translation.Y;
f  = CameraParameters.Intrinsic.FocalLength;

% Pixel size
% assume square pixels.
sensor_width  = CameraParameters.Intrinsic.Sensor.Width;
sensor_height = CameraParameters.Intrinsic.Sensor.Height;

xw_min = -1/f *  (-sensor_width/2 + f * tx/L) * (Z - L) + tx;
xw_max =  -1/f  * (sensor_width/2 + f * tx/L)  * (Z - L) + tx;

yw_min = -1/f *  (-sensor_height/2 + f * ty/L) * (Z - L) + ty;
yw_max =  -1/f *  (sensor_height/2 + f * ty/L) * (Z - L) + ty;

x_world = [xw_min, xw_max];
y_world = [yw_min, yw_max];

end