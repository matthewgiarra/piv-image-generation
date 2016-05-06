function [I1,I2]=microPIVimg(U, V);

%physical parameters (um,us)
%%%%%%%%%%%%%%%%

% Illumination wavelength in microns
lamda = 0.532;

% Objective magnification
M = 40;

% Numerical aperture of objective lens
NA  = 0.4;

% Pixel size in microns
pixel_size_microns  = 19;

% Max intensity
Imax  = 1;

% Diameter of cylinrical channel in microns?
Dc    = 340;

% Challen length in microns?
L     = 340;

% Average velocity in m/s
Vbulk = 0.005773;       %m/s

% Time separation between frames, 
% probably in milliseconds? 
dt    = 1000;

% Particle diameter in microns
dp_microns    = 1; % Particle diameter microns

% Particle Concentration (volume percent?) 
C = 0.005;

% Temperature (Kelvin)
T = 0;

% Fluid viscosity in Pa*S
% Water is 1 cP = 0.001 Pa*S = 1 mPa * S
visc = 1.12E-3;

% Boltzman's constant
% m^2 * kg / (s^2 * K)
kb = 1.38064852E-23; 

% Particle diameter in meters
dp_meters = dp_microns / 10^6;

% Calculate diffusion (m^2 / s)
% The factor of 10^12 converts meters^2 to microns^2
% In the original code this was 10^6 not 10^12, which
% seems wrong to me. 
Diff = kb * T / (6 * pi * visc * dp_meters / 2) * 10^12;

% % Displacements in pixels
% U = 3.00 + 1 * 0.36436;
% V = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%% 

%conversion to pixel units

% Diameter of the channel in pixels
Dc    = round(Dc*M/pixel_size_microns);  % (microns / (microns/pix))) 

% Length of the channel in pix
L     = round(L*M/pixel_size_microns);   

% Factor of 2 gives the max displacement in pix 
Vmax  = (Vbulk*2)*M/pixel_size_microns*dt; 

% Particle diameter in pix 
dp    = dp_microns * M / pixel_size_microns;           

% Diffraction limited diameter in pixels
ds    = 2.44 * (M+1)*lamda / (2 * NA) /pixel_size_microns;  

% Diameter of the particle image in pixels
d     = sqrt(dp^2+ds^2);

% I think this is the f#? because NA = 1 / (2 * f#)
zf    = 1 / (2*NA);

% Number of particles per cubic pixel
C     = C / (pi/6*(dp)^3);

% Laser wavelength in units of pixels
lamda = lamda * M / pixel_size_microns;

% Brownian motion in units of pixels per frame
% This is the standard deviation of the Guassian
% distribution from which random displacements are drawn
BM    = sqrt(2 * Diff * dt) * M / pixel_size_microns;

% % Find number of particles and nearfield depth

% I think this iz the Z position corresponding
% to the limit of the near-field, but 
% I can't figure out why this equation is correct.
zmax = (6 / (pi * C))^(1/3) / NA; 

% I think this is the diameter of a particle
% at the limit of the near field?
Lplus = sqrt(d^2 + (2*NA*zmax)^2);

N = round(C*(L+Lplus+Vmax)*Dc^2);  %image volume + Vmax + Lplus
B = sqrt(8);

fprintf('\n----------------------------------\n');
fprintf('channel diameter    = %0.0f pixels\n',Dc)
fprintf('channel length      = %0.0f pixels\n\n',L)

fprintf('magnification       = %0.0f pixels\n',M)
fprintf('sensor size         = %0.0f pixels\n',pixel_size_microns)
fprintf('um/pix              = %0.2f pixels\n\n',pixel_size_microns/M)

fprintf('number of particles = %0.0f\n\n',N)

fprintf('volume fraction %%   = %0.3f\n',C*(pi/6*dp^3)*100)
fprintf('particle density    = %0.8e\n\n',C)

fprintf('particle diameter   = %0.2f pixels\n',dp)
fprintf('PSF diameter        = %0.2f pixels\n',ds)
fprintf('effective diameter  = %0.2f pixels\n\n',d)

fprintf('nearfield depth     = %0.0f pixels\n',zmax)
fprintf('maximum diameter    = %0.2f pixels\n',Lplus)
fprintf('separation length   = %0.2f pixels\n\n',(3/4/pi/C)^(1/3))

fprintf('bulk velocity       = %0.2f m/s\n',Vbulk)
fprintf('pulse separation    = %0.0f us\n',dt)
fprintf('bulk image shift    = %0.1f pixels\n',Vmax/2)
fprintf('max image shift     = %0.1f pixels\n',Vmax)

fprintf('diffusion           = %0.3f um^2/us\n',Diff)
fprintf('Brownian motion     = %0.3f pix\n',BM)

fprintf('----------------------------------\n');

fprintf('\nGenerating Random Variables...')

%random variable generation
X = (L+Lplus+Vmax)*rand(N,1)-(Lplus/2+Vmax);
Y = Dc*rand(N,1);
Z = Dc*rand(N,1)-Dc/2;

de = sqrt(d^2+(Z./zf).^2);
Io = (de.^-2).*Imax*d^2;

fprintf('done\n')
fprintf('Generating Images...\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Image 1                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve for nearfield particles
I1 = zeros(L+1,Dc+1);
% PP2 = zeros(L+1,Dc+1,N);
% PP3 = zeros(L+1,Dc+1,N);
cp=0;

for p=1:N
    
Particle_Clock = clock;
    
%     if rem(p,500)==0
%         fprintf('%0.0f\n',p)
%     end

    % Only render "in focus" particles
    if abs(Z(p)) < (zmax)
        % Particle counter
        cp=cp+1;
        %fprintf('%0.0f\n',cp)
        
        for i=floor(X(p)-0.75*de(p)):ceil(X(p)+0.75*de(p))
            for j=floor(Y(p)-0.75*de(p)):ceil(Y(p)+0.75*de(p))

                if i>=0 && i<=L && j>=0 && j<=Dc && sqrt((i-X(p))^2+(j-Y(p))^2)<0.75*de(p)

                    %gaussian function value
                    I1(i+1,j+1) = I1(i+1,j+1) + Io(p)*exp(-8*((i-X(p))^2+(j-Y(p))^2)/de(p)^2);

                end

            end
        end
        
    end

% Particle_Elapse = etime(clock,Particle_Clock);
% fprintf('%4.0f of %0.0f \t %0.0f:%5.2f\n',p,N,floor(Particle_Elapse/60),rem(Particle_Elapse,60));    

end

%solve for background glow
y = 0:Dc;
Imask = zeros(L+1,Dc+1);

% What the hell is q1?
for q1=zmax:Dc/2

    
    dez = sqrt(d^2+(q1./zf).^2);
    Deff = Dc;
    
    % What the hell is "Line"?
    Line = pi/4*2*C*Imax*d^2/4*( erf(sqrt(8)*(y-(Dc-Deff)/2)/dez) - erf(sqrt(8)*(y-Dc+(Dc-Deff)/2)/dez) );

        
    for q2=0:L
        Imask(q2+1,:)=Imask(q2+1,:)+Line;
    end

end

I1=I1' + Imask';
fprintf('Image 1 complete\n');

% figure(1),imagesc(I1),axis image,colormap gray,pause(0.01)
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Image 2                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%displace the particles
% U = rectflow(Y-Dc/2,Z,Dc,Dc);

% X = X + Vmax/2*2.1174*U/max(U(:));
X = X + U;
Y = Y + V;
% X = X + shear*Z;

%brownian motion
X = X + BM*randn(N,1);
Y = Y + BM*randn(N,1);
Z = Z + BM*randn(N,1);

%solve for the nearfield particles
I2 = zeros(L+1,Dc+1);
cp=0;

for p=1:N

    if rem(p,500)==0
        fprintf('%0.0f\n',p)
    end

    if abs(Z(p))<zmax

        cp=cp+1;
        %fprintf('%0.0f\n',cp)
        
        for i=floor(X(p)-0.75*de(p)):ceil(X(p)+0.75*de(p))
            for j=floor(Y(p)-0.75*de(p)):ceil(Y(p)+0.75*de(p))

                if i>=0 && i<=L && j>=0 && j<=Dc && sqrt((i-X(p))^2+(j-Y(p))^2)<0.75*de(p)

                    %gaussian function value
                    I2(i+1,j+1) = I2(i+1,j+1) + Io(p)*exp(-8*((i-X(p))^2+(j-Y(p))^2)/de(p)^2);

                end

            end
        end
        
    end
    
end

%solve for the background glow
y = 0 : Dc;
Imask = zeros(L+1,Dc+1);

for q1=zmax:Dc/2

    %fprintf('%0.0f\n',q1)
    dez = sqrt(d^2+(q1./zf).^2);
    Deff = Dc;
    
    % What the hell is "Line"? 
    Line = pi/4*2*C*Imax*d^2/4*( erf(sqrt(8)*(y-(Dc-Deff)/2)/dez) - erf(sqrt(8)*(y-Dc+(Dc-Deff)/2)/dez) );

    % What the hell is q2?
    for q2=0:L
        
        % What the hell is Imask?
        Imask(q2+1,:)=Imask(q2+1,:)+Line;
    end

end

% What the hell is this?
I2 = I2' + Imask';

fprintf('Image 2 complete\n');
% keyboard
% figure(2),imagesc(I2),axis image,colormap gray,pause(0.01)

%Normalized Intensity
% Iimg = max([max(max(I1)) max(max(I2))]);
I1 = I1/2;
I2 = I2/2;

end
