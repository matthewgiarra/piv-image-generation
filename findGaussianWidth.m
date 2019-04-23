function [STDY, STDX] = findGaussianWidth(IMAGESIZEY, IMAGESIZEX, WINDOWSIZEY, WINDOWSIZEX)
% FINDGAUSSIANWIDTH determines the standard deviation of a normalized Gaussian function whose
% area is approximately equal to that of a top-hat function of the desired
% effective window resolution. 
% 
% John says that the volume under a this Gaussian window will not be
% exactly the volume under the square window. Check into this.
% 
% INPUTS
%      xregion = Width of each interrogation region (pixels)
%      yregion = Height of each interrogation region (pixels)
%      xwin = Effective window resolution in the x-direction (pixels)
%      ywin = Effective window resolution in the y-direction (pixels)
% 
% OUTPUTS
%       sx = standard deviation of a normalized Gaussian for the x-dimension of the window (pixels)
%       sy = standard deviation of a normalized Gaussian for the y-dimension of the window (pixels)
% 
% EXAMPLE
%       xregion = 32;
%       yregion = 32;
%       xwin = 16;
%       ywin = 16;
%       [sx sy] = findGaussianWidth(xregion, yregion, xwin, ywin);
% 
% SEE ALSO
% 

% Initial guess for standard deviaitons are half the respective window sizes
% sx = xwin/2;
% sy = ywin/2;

STDX = 50 * WINDOWSIZEX;
STDY = 50 * WINDOWSIZEY;

% Generate x and y domains
x = - IMAGESIZEX/2 : IMAGESIZEX/2;
y = - IMAGESIZEY/2 : IMAGESIZEY/2;

% Generate normalized zero-mean Gaussian windows 
xgauss = exp(-(x).^2/(2 * STDX^2));
ygauss = exp(-(y).^2/(2 * STDY^2));

% Calculate areas under gaussian curves
xarea = trapz(x, xgauss);
yarea = trapz(y, ygauss);

if WINDOWSIZEX < xarea

    % Calculate initial errors of Gaussian windows with respect to desired
    % effective window resolution
    xerr = abs(1 - xarea / WINDOWSIZEX);

    % Initialize max and min values of standard deviation for Gaussian x-window
    sxmax = 100 * IMAGESIZEX;
    sxmin = 0;

    % Iteratively determine the standard deviation that gives the desired
    % effective Gaussian x-window resolution

    % Loop while the error of area under curve is above the specified error tolerance
    while xerr > 1E-5

        %  If the area under the Gaussian curve is less than that of the top-hat window
        if xarea < WINDOWSIZEX
            % Increase the lower bound on the standard deviation
            sxmin = sxmin + (sxmax - sxmin) / 2;
        else
            % Otherwise, increase the upper bound on the standard deviation
            sxmax =  sxmin + (sxmax - sxmin) / 2;
        end

        % Set the standard deviation to halfway between its lower and upper bounds
        STDX = sxmin + (sxmax - sxmin) / 2;

        % Generate a Gaussian curve with the specified standard deviation
        xgauss = exp(-(x).^2/(2 * STDX^2));

        % Calculate the area under this Gaussian curve via numerical
        % integration using the Trapezoidal rule
        xarea = trapz(x,xgauss);

        % Calculate the error of the area under the Gaussian curve with
        % respect to the desired area 
        xerr = abs(1 - xarea / WINDOWSIZEX);

    end

end

if WINDOWSIZEY < yarea
    
    yerr = abs(1 - yarea / WINDOWSIZEY);

    % Initialize max and min values of standard deviation for Gaussian y-window
    symax = 100 * IMAGESIZEY;
    symin = 0;

    % Iteratively determine the standard deviation that gives the desired
    % effective Gaussian y-window resolution

    % Loop while the error of area under curve is above the specified error tolerance
    while yerr > 1E-5

        % If the area under the Gaussian curve is less than that of the top-hat window
        if yarea < WINDOWSIZEY
            % Increase the lower bound on the standard deviation
            symin = symin + (symax - symin) / 2;
        else
            % Otherwise, increase the upper bound on the standard deviation
            symax =  symin + (symax - symin) / 2;
        end

        % Set the standard deviation to halfway between its lower and upper bounds
        STDY = symin + (symax - symin) / 2;

        % Generate a Gaussian curve with the specified standard deviation
        ygauss = exp(-(y).^2/(2 * STDY^2));

        % Calculate the area under this Gaussian curve via 
        % numerical integration using the Trapezoidal rule
        yarea = trapz(y, ygauss);

        % Calculate the error of the area under the Gaussian curve with
        % respect to the desired area 
        yerr = abs(1 - yarea / WINDOWSIZEY);

    end

end

end








