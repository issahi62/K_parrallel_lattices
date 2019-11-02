% INITIALIZE CODES 
clear all
close all
clc


%% Shaping a matrix 

% PARAMETERS 

Sy = 3e6; % physical size 
Sx = Sy; 

Nx = 100; % Number of cells in x
Ny = Nx;

%% Define width of rectangle 
wx = 0.2; 
wy = 0.6; 

% Compute position indices 
dx = Sx/Nx; 
dy = Sy/Ny; 

xa = [0:Nx-1]*dx; xa = 2*(xa-mean(xa)); 
ya = [0:Ny-1]*dy; ya = 2*(ya-mean(ya));

kx = xa; 
ky = ya; 

[Y, X] = meshgrid(ky,kx); 

% DEFINE TWO POINTS ALONG BAR AND WIDTH
w =0.1;
Z=1;
r1 = [ 0.1 ; 0 ; 0.8 ]; 
r2 = [ 0.8 ; 1 ; 0.1 ];
% CALCULATE DISTANCE FUNCTION
d =r2-r1;
d = d/norm(d);
D = sqrt( (d(2)*(r1(3) - Z) - d(3)*(r1(2) - Y)).^2 + ...
(d(3)*(r1(1) - X) - d(1)*(r1(3) - Z)).^2 + ... 
(d(1)*(r1(2) - Y) - d(2)*(r1(1) - X)).^2 );
% CALCULATE BAR
ER = (D <= w);

% DEFINE TWO POINTS
x1 = -0.50; y1 = +0.25; x2 = +0.50; y2 = -0.25;
% FILL HALF SPACE
m = (y2 - y1)/(x2 - x1);
A = (Y - y1) - m*(X - x1) > 0;

 % program calculates temperature distribution
clear % clears variables
clear global % clears global variables
format short
Nx = 3;% number of points along x axis
Ny = 3;% number of points along y axis
N = Nx*Ny;%total number of points
for i = 1:N% setting off diagonals
a1(i) = 1;
a2(i) = 1;
a4(i) = 1;
a5(i) = 1; end
for i = 1:Ny% setting zeros in off diagonals
a4((i-1)*Nx+1) = 0;
a2(i*Nx) = 0;
end
for i = 1:N% setting the main diagonal
a3(i) = -4;
end
for i = 1:N% setting the augmented vector
b(i) = -300;
end
b(1) = -600;
b(Nx) = -600;
b(N-Nx+1) = -600;
b(N) = -600;
for j = 2:Ny-1
for i = 2:Nx-1
       b((j-1)*Nx+i) = 0;
end
end
c = [-Nx,-1,0,1,Nx];%setting the matrix pattern
%forming directly the matrix in sparse format
A = spdiags([a1' a2' a3' a4' a5'], c, N);
% a1-a5 are transposed
x = A\b;% solving the set of equations







