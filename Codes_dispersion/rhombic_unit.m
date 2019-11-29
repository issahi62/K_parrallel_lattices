%*******************
%% UNIT CELL
%****************
L1=0.3; 
L2 =0.4;
a=5;
Nx=512;
Ny=512;
Mx=20;               %% map X [m]
My=20;               %% map Y [m]

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
idxa = X.^2+ Y.^2 >.01;


x=linspace(-Mx/2,Mx/2,Nx);
y=linspace(-My/2,My/2,Ny);
[X,Y]=meshgrid(x,y);
idx1=  (abs(X)<a*sqrt(3)/2);
idx2=(tan(pi/6)*X+a>Y) .* (tan(pi/6)*X-a<Y);
idx=idx1.*idx2;
ER = idx.*idxa; 

%*******************
%% PLOT
%****************
figure(4); 
pcolor(ER)
shading interp; 
%colormap(jet);