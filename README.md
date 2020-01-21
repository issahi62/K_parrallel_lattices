##### SOURCE BEGIN #####
% PWEM developed by Ibrahim Issah.

% INITIALIZE MATLAB 
 close all; 
 clc; 
 clear; 
 
% 1 = square_lattice 
% 2 = rectangular_lattice
% 3 = rhombic_lattice
% 4 = DFB_unit
% 5 = oblique shape
% 6 = ellipsoid 



type = 1; 


% OPEN FIGURE WINDOW 
figure('Color', 'w', 'Position', [360 278 807 420]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
field_selector = 0; 
% UNIT -CELL 

a = 1;
r = 0.15*a; %radius of cylinder
er = 9.0; %permittivity

% MODE SELECTOR 
mm = 2;
mm_selector = field_selector; 

% HIGH RESOLUTION GRID 
Nx = 512; 
Ny = Nx; 

% PWEM PARAMETERS 
P=7;  % ensure this is an odd number 
Q = P; 

NBx = 20; 
NBy= NBx; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD UNIT CELL
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MESHGRID TO Form the Unit cell
dx = a/Nx; 
dy = a/Ny; 

xa = [1:Nx]*dx; xa = xa-mean(xa); 
ya = [1:Ny]*dy; ya = ya-mean(ya); 

[Y,X] = meshgrid(ya,xa);

% BUILD UNIT CELL
UR = ones(Nx,Ny); % permittivity;
%********************************
%% SELECT A UNIT CELL OF YOUR CHOICE
%%
%*******************************
switch (type)
    case 1
        square_unit;
    case 2
        rectangular_unit;
    case 3
        rhombic_unit;
    case 4
        DFB_unit;
    case 5
        oblique;
    otherwise
        ellip_unit;
end 


ER = 1 + (er-1.0)*ER; 
%ER = ones(Nx,Ny); 
% COMPUTE CONVOLUTION MATRICES 
URC = convmat_PWMEM(UR, P, Q); 
ERC = convmat_PWMEM(ER, P, Q); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE MESHGRID OF BLOCH WAVE VECTORS
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betax = (pi/a)*(linspace(-1, 1, NBx)); 
betay = (pi/a)*(linspace(-1, 1, NBy)); 

[BY, BX] = meshgrid(betay, betax);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM PWEM
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE SPATIAL HARMONIC INDICES 
p= [-floor(P/2): floor(P/2)]; %indices along x
q= [-floor(Q/2): floor(Q/2)]; % indices along y

% INITIALIZE BAND DATA 
KO = NaN*ones(NBx, NBy, P*Q); 
KOH = NaN*ones(NBx, NBy, P*Q);
Emode= NaN*ones(NBx, NBy, P*Q);
Hmode= NaN*ones(NBx, NBy, P*Q);
Vsto = NaN*ones(2*NBx+9, 2*NBy+9, NBy);
%% INSERTION OF FOURIER COEFFICEINTS INTO LARGE GRID
%%
nxc = ceil(Nx/2); 
nx1 = nxc-floor(P/2); 
nx2 = nxc + floor(P/2); 
nyc= ceil(Ny/2); 
ny1 = nyc-floor(P/2); 
ny2 = nyc + floor(P/2); 

sf = zeros(Nx,Ny); 
%
% MAIN LOOP REPLACE_WITH_DASH_DASH ITERATE OVER BLOCH WAVE VECTOR
%

for nby = 1:NBy
    for nbx = 1:NBy
    
        %GET NEX BETA 
        bx = BX(nbx, nby);
        by = BY(nbx, nby); 

        % FORM K MATRICES
%%
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
        %Emode(nbx, nby, :) = sort(diag(fftshift(fft(V)))); 
        
         Vsto(:, :, nbx) = V; 

        % COMPUTE EIGEN-VALUES FOR H MODE 
        A = KX/ERC*KX + KY/ERC*KY; 
        B = URC; 

        [VH,DH] = eig(full(A), full(B)); 
        DH = real(sqrt(diag(DH)));
        KOH(nbx, nby, :) = sort(DH);  
        %Hmode(nbx, nby, :) = sort(diag(fftshift(fft(VH)))); 
    end   
     %drawnow; 
end

% CALCULATE NORMALIZED FREQUENCY 
WN = a*KO/(2*pi);
WNH = a*KOH/(2*pi);
%[px,py] = gradient(WN(:, :, mm));
%% VISUALIZATION OF FIELDS
%%
%s = VH(:, m); 
  storer = zeros(Nx, Ny, P*Q);
%% COMPUTE COMPLEX E-FIELD FROM EIGEN VECTORS.
%%
        for kobby = 1:(P*Q)
        s = Vsto(:,kobby,2);
        s = reshape(s, P, Q); 
        sf(nx1:nx2, ny1:ny2) =s;
        az = ifft2(ifftshift(sf)); 
        az = az/max(abs(az(:))); 
        phase = exp(-1i*(betax(1)*X + betay(2)*Y)); 
        Ez = phase.*az;
        storer(:, :, kobby) =Ez; 
        end
%[px1,py1] = gradient( storer(:, :, mm),-1, -1);
% E MODE
 figure(1); 
    for m = 1:9 % showing only 6 bands 
        subplot(3,3,m); 
        pcolor((WN(:, :, m+mm_selector).')); hold on
        %quiver(betax,betay,px,py), hold off, axis image
        shading interp
        axis equal tight 
        colorbar 
        colormap(jet);
        title(['E Mode: ' num2str(m+mm_selector)]); 
    end
figure('Color', 'w', 'Position', [360 278 807 420]);
figure(2); 
    for m = 1:9 % showing only 6 bands 
        subplot(3,3,m);
        
        pcolor(betax, betay.', WNH(:, :, m+mm_selector).'); hold on 
        %quiver(betax,betay, px, py), hold off, axis image
        shading interp
        axis equal tight 
        colorbar 
        colormap(jet);
        title(['H Mode: ' num2str(m+mm_selector)]); 
    end
    
    figure('Color', 'w', 'Position', [360 278 807 420]);
    figure(5); 
    
    for m = 1:9 % showing only 6 bands 
        subplot(3,3,m);
        
        pcolor(xa, ya.', abs(storer(:, :, m+mm_selector).').^2); hold on 
       % quiver(xa,ya.',real(px1*cos(45)),real(py1*sin(45))), hold off, axis image
        shading interp
        axis equal tight 
        colorbar 
        colormap(jet);
        title(['E Field: ' num2str(m+mm_selector)]); 
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DRAW ISO-FREQUENCY CONTOURS
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLEAR FIGURE WINDOW 
%clf;

% CONTOURS
figure(3); 
[c, h]= contour(betax, betay, WN(:, :, mm), 'LineWidth', 2.5); hold on; 
%quiver(betax,betay,px,py), hold off, axis image
clabel(c, h); 
title(['Iso-frequency of mode: ', num2str(mm)]); 
colormap(hot)

% SET VIEWS
axis equal tight 
xlabel('$\beta_x$', 'Interpreter', 'Latex'); 
ylabel('$\beta_y$', 'Interpreter', 'Latex', 'Rotation', 0, ...
    'HorizontalAlignment', 'right');
##### SOURCE END #####
--></body></html>
