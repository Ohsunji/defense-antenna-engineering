function [ E_x, E_y, E_z, x, y ] = BackProjection_PlanarNearField( E_meas_x, E_meas_y, z_0, f, dx, dy )
% function [ E_x, E_y, E_z, x, y ] = BackProjection_PlanarNearField( E_meas_x, E_meas_y, z_0, f, dx, dy )
% AUTHOR: B*log10(Engineering)
% WEBSITE: http://blog.p86.dk/
% DATE: 25 Jan 2012
%
% DESCRIPTION:
% This function can used for antenna diagnostics, by taking planar near-field
% measured data, i.e. the vertically and horizonally measured tangential 
% E-fields, and back projecting them to the E-fields across the aperture 
% of the antenna under test (AUT). To prevent aliasing spatial sampling
% steps should be less than half the wave length.
%
% Note that often an open ended waveguide is used as a probe for planar 
% measurements, therefore it is expected that probe correction has been 
% applied to the measured data before using this function.
%
% This script is based on the 5 steps described by Balanis in [1, pg. 706].
% Inspiration for writing this script was also taken from [2].
% 
% PARMAETERS:
%   E_meas_x:   Planar Near-Field Horizontally Measured Data
%   E_meas_y:   Planar Near-Field Vertically Measured Data
%   z_0:        Distance between probe and AUT [m]
%   f:          Frequency [Hz]
%   dx:         Planar x-plane step size [m]
%   dy:         Planar y-plane step size [m]
%
% RETURNS:
%   E_x:        Back projected E_x field
%   E_y:        Back projected E_y field
%   E_z:        Back projected E_z field
%   x:          Planar x-axis [m] (used for plotting)
%   y:          Planar y-axis [m] (used for plotting)
%
% REFERENCES:
%   [1] Balanis C.A., 2005, Antenna Theory Analysis and Design, 3rd Edition, 
%       John Wiley & Sons, Chapter 12.9
%   [2] Van Caekenberghe K., Logan J., Mynster A.P., Pelk M.J., Ponder C.,
%       http://www.mathworks.com/matlabcentral/fileexchange/23385-nf2ff/content/NF2FF.m,
%       Viewed: 24 Jan 2012

% Intialise Variables
E_x = 0;
E_y = 0;
E_z = 0;
x = 0;
y = 0;

% Check X and Y Measured Data Matrices are the same size
if size(E_x) ~= size(E_y)
   disp('Error! In "BackProjection_PlanarNearField": E_x and E_y are not the same size.');
   return;
end

c = 299792458; % Speed of light in vacuum [m/s]
lambda0 = c/f;   % Freespace Wavelength [m]
k0 = 2*pi/lambda0; % Freespace Wavenumber

% Check if step sizes will introduce aliasing
if lambda0/2 < dx
   disp('Warning! X-Axis step size in larger than half the wave length, which will lead to aliasing.');
end

if lambda0/2 < dy
    disp('Warning! Y-Axis step size in larger than half the wave length, which will lead to aliasing.');
end

S = size(E_meas_x);
M = S(1);   % Number of x-axis samples
N = S(2);   % Number of y-axis samples

a = dx*(M-1); % Length of the Scan Plane (x direction) [m]
b = dy*(N-1); % Height of Scan Plane (y direction) [m]
x = -a/2:a/(M-1):a/2; % Scan plane x-axis points
y = -b/2:b/(N-1):b/2; % Scan plane y-axis points

MI = 10*M;
NI = 10*N;
m = (-MI/2):1:(MI/2-1);
n = (-NI/2):1:(NI/2-1);

kx = 2*pi*m/(MI*dx);   % Define kx axis for plane wave spectrum
ky = 2*pi*n/(NI*dy);   % Define ky axis for plane wave spectrum
[ky_grid, kx_grid] = meshgrid(ky,kx);
kz_grid = sqrt(k0^2 - kx_grid.^2 - ky_grid.^2); % Define kz axis for plane wave spectrum
  
% 2D Fourier Transform Measured Data to Plane Wave Spectrum
% Note: A different time convention is used therefore ifft2 is used 
% instead of fft2.
fx = ifftshift(ifft2(E_meas_x,MI,NI)); % See Eq. (12-85a)
fy = ifftshift(ifft2(E_meas_y,MI,NI)); % See Eq. (12-85b)
% kz axis could also be define but can be found from kx and ky, see Eq. (12-84b)

% Perform Back Projection, by z_0
fx_z0 = zeros(MI,NI); % Allocate Memory
fy_z0 = zeros(MI,NI); % Allocate Memory
fz_z0 = zeros(MI,NI); % Allocate Memory
for iy = 1:1:NI
    for ix = 1:1:MI
        % Propogating Waves
        if(isreal(kz_grid(ix,iy)))
            fx_z0(ix,iy) = fx(ix,iy)*exp(1i*kz_grid(ix,iy)*z_0); % See Eq. (12-86)
            fy_z0(ix,iy) = fy(ix,iy)*exp(1i*kz_grid(ix,iy)*z_0); % See Eq. (12-86)
            fz_z0(ix,iy) = -(fx_z0(ix,iy)*kx_grid(ix,iy) + fy_z0(ix,iy)*ky_grid(ix,iy))/kz_grid(ix,iy); % See Eq. (12-86a)
            
        % Evanescent Waves
        else
            fx_z0(ix,iy) = 0;
            fy_z0(ix,iy) = 0;
            fz_z0(ix,iy) = 0;
        end
    end
end

% Perform Inverse Fourier Transform of the Back Projected Plane Wave Spectrum
% Note: A different time convention is used therefore fft2 is used 
% instead of ifft2.
E_x_Zero_Padded = fft2(ifftshift(fx_z0)); % See Eq. (12-86)
E_y_Zero_Padded = fft2(ifftshift(fy_z0)); % See Eq. (12-86)
E_z_Zero_Padded = fft2(ifftshift(fz_z0)); % See Eq. (12-86)

% Crop Data
E_x = E_x_Zero_Padded(1:M,1:N);
E_y = E_y_Zero_Padded(1:M,1:N);
E_z = E_z_Zero_Padded(1:M,1:N);
