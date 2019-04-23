function WINDOW = gaussianWindowFilter(DIMENSIONS, WINDOWSIZE, WINDOWTYPE)
% gaussianWindowFilter(DIMENSIONS, WINDOWSIZE, WINDOWTYPE) creates a 2-D
% gaussian window
%
% INPUTS
%   DIMENSIONS = 2 x 1 Vector specifying  the Dimensions (in rows and columns) of the gaussian window filter
%   WINDOWSIZE = 2 x 1 vector specifying the effective window resolution of
%   the gaussian window. WINDOWSIZE can either specify the resolution in
%   pixels or as a fraction of the filter dimensions. This option is
%   controlled by the input WINDOWTYPE.
%   WINDOWTYPE = String specifying whether WINDOWSIZE specifies a
%   resolution in pixels ('pixels') or as a fraction of the window
%   dimensions ('fraction').
% 
% OUTPUTS
%   WINDOW = 2-D gaussian window
% 
% SEE ALSO
%   findwidth

% Default to an absolute size window type
if nargin < 3
    WINDOWTYPE = 'fraction';
end

if length(DIMENSIONS) == 1
    dims = DIMENSIONS * [1, 1];
else
    dims = DIMENSIONS;
end

% Signal height and width
height = dims(1);
width = dims(2);

if length(WINDOWSIZE) == 1
    win_size = WINDOWSIZE * [1, 1];
else
    win_size = WINDOWSIZE;
end

% Determine whether window size is an absolute size or a fraction of the
% window dimensions
if strcmp(WINDOWTYPE, 'fraction')
    windowSizeX = width .* win_size(2);
    windowSizeY = height .* win_size(1);
elseif strcmp(WINDOWTYPE, 'pixels');
    windowSizeX = win_size(2);
    windowSizeY = win_size(1);
else
    error('Invalid window type "%s"\n', WINDOWTYPE);
end

% Standard deviations
[sy, sx] = findGaussianWidth(height, width, windowSizeY, windowSizeX);

% Calculate center of signal
xc = (width-1)/2;
yc = (height-1)/2;

% Create grid of x,y positions to hold gaussian filter data
[xo,yo] = meshgrid(0:(width-1), 0:(height-1));

% Shift the coordinates to make them 
% symmetric about the centroid of the array
x = xo - xc;
y = yo - yc;

% Calculate gaussian distribution (X)
WindowX = exp( - (x.^2 / (2 * sx^2)));

% Calculate gaussian distribution (Y)
WindowY = exp( - (y.^2 / (2 * sy^2)));

% 2-D Gaussian Distribution
WINDOW = WindowX .* WindowY;


end