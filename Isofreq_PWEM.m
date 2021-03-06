% demo_bandsPWEM

% INITIALIZE MATLAB 
close all; 
clc; 
clear all; 

% OPEN FIGURE WINDOW 
figure('Color', 'w', 'Position', [360 278 807 420]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT -CELL 

a = 1;
r = 0.39*a; %radius of cylinder
er = 9.0; %permittivity

% MODE SELECTOR 
mm = 2; 

% HIGH RESOLUTION GRID 
Nx = 512; 
Ny = Nx; 

% PWEM PARAMETERS 
P=11;  % ensure this is an odd number 
Q = P; 

NBx = 20; 
NBy= NBx; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD UNIT CELL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MESHGRID TO Form the Unit cell
dx = a/Nx; 
dy = a/Ny; 

xa = [1:Nx]*dx; xa = xa-mean(xa); 
ya = [1:Ny]*dy; ya = ya-mean(ya); 

[Y,X] = meshgrid(ya,xa);

% BUILD UNIT CELL
UR = ones(Nx,Ny); % permittivity;
ER = (X.^2+Y.^2)>r^2; 
ER = 1.0 + (er-1.0)*ER; 

% COMPUTE CONVOLUTION MATRICES 
URC = convmat_PWMEM(UR, P, Q); 
ERC = convmat_PWMEM(ER, P, Q); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE MESHGRID OF BLOCH WAVE VECTORS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betax = (pi/a)*(linspace(-1, 1, NBx)); 
betay = (pi/a)*(linspace(-1, 1, NBy)); 

[BY, BX] = meshgrid(betay, betax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM PWEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE SPATIAL HARMONIC INDICES 
p= [-floor(P/2): floor(P/2)]; %indices along x
q= [-floor(Q/2): floor(Q/2)]; % indices along y

% INITIALIZE BAND DATA 
KO = NaN*ones(NBx, NBy, P*Q); 
Emode= NaN*ones(NBx, NBy, P*Q);

%
% MAIN LOOP -- ITERATE OVER BLOCH WAVE VECTOR
%

for nby = 1:NBy
    for nbx = 1:NBy
    
        %GET NEX BETA 
        bx = BX(nbx, nby);
        by = BY(nbx, nby); 

        % FORM K MATRICES
        KX = bx-2*pi*p/a; 
        KY = by - 2*pi*q/a; 
        [KY, KX] = meshgrid(KY, KX); 
        KX = diag(sparse(KX(:))); 
        KY = diag(sparse(KY(:))); 

        % COMPUTE EIGEN-VALUES FOR E MODE 
        A = KX/URC*KX + KY/URC*KY; 
        B = ERC; 

        [V,D] = eig(full(A), full(B)); 
        D = real(sqrt(diag(D)));

        KO(nbx, nby, :) = sort(D);  
        Emode(nbx, nby, :) = sort(diag(V)); 
    end 
    % SHOW BANDS 
    figure(1); 
    for m = 1:9 % showing only 6 bands 
        subplot(3,3,m); 
        pcolor(betax, betay.', KO(:, :, m).')
        shading interp
        axis equal tight 
        colorbar 
        colormap(hot);
        title(['Mode: ' num2str(m)]); 
    end
     drawnow; 
end

% CALCULATE NORMALIZED FREQUENCY 
WN = a*KO/(2*pi);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DRAW ISO-FREQUENCY CONTOURS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR FIGURE WINDOW 
%clf;

% CONTOURS
figure(2); 
[c, h]= contour(betax, betay, WN(:, :, mm), 'LineWidth', 2.5);
clabel(c, h); 
title(['Iso-frequency of mode: ', num2str(mm)]); 
colormap(hot)

% SET VIEWS
axis equal tight 
xlabel('$\beta_x$', 'Interpreter', 'Latex'); 
ylabel('$\beta_y$', 'Interpreter', 'Latex', 'Rotation', 0, ...
    'HorizontalAlignment', 'right'); 

% E_mode 
figure(3);
pcolor(abs(Emode(:, :, mm).^2))
shading interp
colormap(hot)
