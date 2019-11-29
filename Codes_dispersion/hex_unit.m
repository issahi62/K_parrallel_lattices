Mx=20e-9;               %% map X [m]
My=20e-9;               %% map Y [m]

x=linspace(-Mx/2,Mx/2,Nx);
y=linspace(-My/2,My/2,Ny);
[X,Y]=meshgrid(x,y);
idx1=  (abs(X)<a*sqrt(3)/2);
idx2=(tan(pi/6)*X+a>Y) .* (tan(pi/6)*X-a<Y) .* (-tan(pi/6)*X-a<Y) .* (-tan(pi/6)*X+a>Y);
idx=idx1.*idx2;
ER = idx.*adx; 

%*******************
%% PLOT
%****************
figure(4); 
pcolor(ya, xa, ER)
shading interp; 
colormap(jet);