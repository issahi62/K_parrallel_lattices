% demo_bandsPWEM

% INITIALIZE MATLAB 
 close all; 
 clc; 
 clear all; 
 
% 1 = square_lattice 
% 2 = rectangular_lattice
% 3 = ellipsoid 

type = 3; 


% OPEN FIGURE WINDOW 
figure('Color', 'w', 'Position', [360 278 807 420]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
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
P=3;  % ensure this is an odd number 
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
%********************************
%% SELECT A UNIT CELL OF YOUR CHOICE 
%*******************************
switch (type)
    case 1
        square_unit;
    case 2
        rectangular_unit;
    otherwise
        ellip_unit;
end 


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
        Emode(nbx, nby, :) = sort(diag(fftshift(fft(V)))); 
    end 
    
     %drawnow; 
end
% SHOW BANDS 
%     figure(1); 
%     for m = 1:9 % showing only 6 bands 
%         subplot(3,3,m); 
%         pcolor(betax, betay.', KO(:, :, m+mm_selector).')
%         shading interp
%         axis equal tight 
%         colorbar 
%         colormap(jet);
%         title(['Mode: ' num2str(m+mm_selector)]); 
%    end
% CALCULATE NORMALIZED FREQUENCY 
WN = a*KO/(2*pi);
 figure(1); 
    for m = 1:9 % showing only 6 bands 
        subplot(3,3,m); 
        pcolor(betax, betay.', WN(:, :, m+mm_selector).')
        shading interp
        axis equal tight 
        colorbar 
        colormap(jet);
        title(['Mode: ' num2str(m+mm_selector)]); 
    end


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
    for m = 1:9 % showing only 6 bands 
        subplot(3,3,m); 
        pcolor(betax, betay.', abs(Emode(:, :, m+field_selector)).')
        shading interp
        axis equal tight 
        colorbar 
        colormap(jet);
        title(['E-Field: ' num2str(m+field_selector)]); 
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EIGEN - VECTOR V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=1:n
%     PSI = reshape(V(:,i),[Ny,Nx]);
%     PSI = invFFT2D(PSI,Ny,Nx)/(dxx*dyy) ;
%     psi_temp = interp2(XX,YY,PSI,X,Y);
%     psi(:,:,i) = psi_temp / sqrt( trapz( y' , trapz(x,abs(psi_temp).^2 ,2) , 1 )  );  % normalisation of the wave function psi
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % here is a small patch due to differences between Octave and Matlab
% % Matlab order the eigen values while Octave reverse it
% 
% if E(1)>E(2)
%   psi=psi(:,:,end:-1:1);
%   E=E(end:-1:1);
% end
% 
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [Vxy] = invFFT2D(Vk2D,Ny,Nx)
% 
% Nkx=length(Vk2D(1,:));
% Nky=length(Vk2D(:,1));
% 
% Nx1=Nx/2-floor(Nkx/2);
% Nx2=Nx/2+ceil(Nkx/2);
% Ny1=Ny/2-floor(Nky/2);
% Ny2=Ny/2+ceil(Nky/2);
% 
% Vk2D00=zeros(Ny,Nx);
% Vk2D00( Ny1+1:Ny2 , Nx1+1:Nx2)=Vk2D;
% Vxy=ifft2(ifftshift(Vk2D00));
% 
% end
