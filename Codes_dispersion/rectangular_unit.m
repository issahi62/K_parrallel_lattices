%*******************
%% UNIT CELL
%****************
L1=0.3; 
L2 =0.4;


%*******************
%% DASHBOARD
%****************
% Nx = 200; 
% Ny = 200; 
dx = 1/Nx; 
dy = 1/Ny; 

xa = [0:Nx-1]*dx; xa = xa-mean(xa); 
ya = [0:Ny-1]*dy; ya = ya-mean(ya);
x0=0;y0=0;     
[X, Y] = meshgrid(ya, xa);

%*******************
%% UNIT CELL
%****************
idx= ( abs(X-x0)<L1 ) .* ( abs(Y-y0)<L2 ) ;
idxa = X.^2+ Y.^2 >r^2;
 
 
ER = idxa.*idx;

%*******************
%% PLOT
%****************
figure(4); 
pcolor(ya, xa, ER)
shading interp; 
colormap('jet'); 