% function G=PSF(row,column, pixel_size, wl, z)  
% %generates the FT of point spread function of NxM size;
% %wavelength lambda at a distance z
% tx=(1:column)-floor(column/2); %define number of columns
% ty=(1:row)-floor(row/2); %define number of rows
% Kx = tx/(pixel_size*column);
% Ky = ty'/(pixel_size*row);
% Kx = ones(row, 1)*Kx;
% Ky = Ky*ones(1, column);
% %% k is outside of sqrt
% % -----
% %k_z=(2*pi/wl)*sqrt(1-(wl*Kx).^2-(wl*Ky).^2);    %z-component of the wave number vector
% %G=exp(1i*k_z*z);%FT of free space point spread function
% G = exp(1i*pi*wl*z*(Kx.^2+Ky.^2));
% % -----


function G=PSF(row,column, pixel_size, wl, n_m,z)  
%generates the FT of point spread function of NxM size;
%wavelength lambda at a distance z
tx=(1-1:column-1)-floor(column/2); %define number of columns
ty=(1-1:row-1)-floor(row/2); %define number of rows
Kx = tx/(pixel_size*column);
Ky = ty'/(pixel_size*row);
Kx = ones(row, 1)*Kx;
Ky = Ky*ones(1, column);
%% k is outside of sqrt
% -----
a = 1-(wl*Kx/n_m).^2-(wl*Ky/n_m).^2;
a = a.*(a > 0);
%k_z=(2*pi/wl)*sqrt(1-(wl*Kx).^2-(wl*Ky).^2);    %z-component of the wave number vector
k_z = (2*pi*n_m/wl)*sqrt(a);
G=exp(1i*k_z*z);%FT of free space point spread function
%G = exp(1i*pi*wl*z*(Kx.^2+Ky.^2));
% -----