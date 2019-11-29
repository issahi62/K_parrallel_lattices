%*******************
%% UNIT CELL
%****************
L1=0.3; 
L2 =0.4;

Rx=4;       %% radius in the the x-direction of the ellipse [m]
Ry=8;       %% radius in the the y-direction of the ellipse [m]
x0=0;y0=0;

%*******************
%% DASHBOARD
%****************
%Nx = 200; 
% Ny = 200; 
dx = 1/Nx; 
dy = 1/Ny; 

xa = [0:Nx-1]*dx; xa = xa-mean(xa); 
ya = [0:Ny-1]*dy; ya = ya-mean(ya);
    
[X, Y] = meshgrid(ya, xa);
     %% center positions of the ellipse [m]
idx= ( abs(X-x0)<L1 ) .* ( abs(Y-y0)<L2 ) ;
idxa=((X-x0)/Rx).^2 + ((Y-y0)/Ry).^2 > .001;


ER = idxa.*idx;

%*******************
%% PLOT
%****************
figure(4); 
pcolor(ya, xa, ER)
shading interp; 