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

a = 100e-9;
r = 0.39*a; %radius of cylinder
er =  12.89; %permittivity 

% HIGH RESOLUTION GRID 
Nx = 512; 
Ny = Nx; 

% PWEM PARAMETERS 
P=31;  % ensure this is an odd number 
Q = P; 
NG2X = 40; 
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

% SHOW UNIT CELL
subplot(1, 3, 1); 
imagesc(xa, ya, ER'); 
colorbar; 
axis equal tight; 

% COMPUTE CONVOLUTION MATRICES 
URC = convmat_PWMEM(UR, P, Q); 
ERC = convmat_PWMEM(ER, P, Q); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE LIST OF BLOCH WAVE VECTORS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RECIPROCAL LATTICE VECTORS
T1 = (2*pi/a) * [1;0]; %  for the reduced brillioun reciprocal space 1
T2 = (2*pi/a) * [0;1]; % y- axis

% KEY POINT OF SYMMETRY 
G = [0;0]; 
X = 0.5*T1; 
M = 0.5*T1 + 0.5*T2; 

% GENERATE LIST 
L1 = norm(G-X); % To find the length from Gamma to X to M in the reduced
L2 = norm(M-G); % Brillioun zone of a square lattice. 
L3 = norm(M-X); 

N1 = NG2X; 
N2 = round(N1*L2/L1); 
N3 = round(N1*L3/L1);

% BLOCH WAVE VECTORS FOR X AXIS AND Y AXIS
BX = [linspace(X(1), G(1), N1), ...
    linspace(G(1), M(1), N2), ...
    linspace(M(1), X(1), N3)]; 

BY = [linspace(X(2), G(2), N1), ...
    linspace(G(2), M(2), N2), ...
    linspace(M(2), X(2), N3)]; 
% 

BETA = [BX ; BY]; 
BETA(:, [N1+1, N1+N2+1]) =[]; 

KP =[ 1, N1, N1+N2-1, N1+N2+N3-2]; % keynotes

KL = { 'X','\Gamma',  'M', 'X'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM PWEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE SPATIAL HARMONIC INDICES 
p= [-floor(P/2): floor(P/2)]; %indices along x
q= [-floor(Q/2): floor(Q/2)]; % indices along y

% INITIALIZE BAND DATA 
NBETA = length(BETA(1,:)); 
KO = zeros(P*Q, NBETA); 
Emode = zeros(P*Q, NBETA);
%
% MAIN LOOP -- ITERATE OVER BLOCH WAVE VECTOR
%

for nbeta = 1:NBETA
    
    %GET NEXT BETA 
    bx = BETA(1, nbeta); 
    by = BETA(2, nbeta); 
    
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
    
    KO(:, nbeta) = sort(D); 
    Emode(:, nbeta) = sort(diag(V));
    
    % SHOW BANDS 
    subplot(1,3,2:3); 
    plot([1:nbeta], KO(:, 1:nbeta), '.b', 'LineWidth', 2); 
    xlim([1 NBETA]); 
    drawnow; 
end 

% CALCULATE NORMALIZED FREQUENCY 
WN = a*KO/(2*pi);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DRAW BAND DIAGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLEAR FIGURE WINDOW 

% PLOT BANDS
plot([1:NBETA], WN, 'LineWidth', 2.5); 

% SET AXIS LIMITs
xlim([1 NBETA]); 
ylim([0, 1]); 

% SET TICK MARKS 
set(gca, 'XTick',KP,'XTickLabel', KL); 

for n = 1: length(KP)
    line(KP(n)*[1 1], [0 1], 'Color', 'k', 'LineStyle', ':')
end 

% LABEL AXES 
xlabel('Bloch Wave VErctor, $\vec{\beta}$', 'Interpreter', 'LaTEx'); 
ylabel('Normalized Frequency, $\omega$', 'Interpreter', 'LaTEx'); 