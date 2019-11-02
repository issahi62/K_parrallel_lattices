%%  Geometry dependence of surface lattice resonances in plasmonic nanoparticle arrays
%% INITIALIZE
close all 
clc
clear all
%******************
%******************
%% DASHBOARD 
%******************
% speed of light 
c = 3.0e8; 
% momentum for periodicity 375nm
px = 415e-9;
py = 375e-9;

% refractive index of the material
n = 1.52; 
% Assumed charge
e1 = 1; 
% planck's constant in elctron_Volt(eV) 
hev = 4.135e-15/(2*pi); 
% Angle of incidence of the incoming light source
theta = pi/2; 
% ky space in the reciprocal space
S1 = 3e6; 
integer_multiple =1; 
N_points = 100; 

%******************
%% OPEN FIGURE 
figure('Color', 'w');
%******************

%******************
%% KY and KX values 
%******************
dx1 = 2*S1/N_points;
ky = [0:N_points-1]*dx1; ky = ky-mean(ky);
kx = ky; 

%******************
%% G parameters 
%******************
G_parameter_x =  integer_multiple*2*pi/px;
G_parameter_y = integer_multiple*2*pi/py;

%******************
%% TE Mode 
%****************** 
E = ((hev*c)/n).*abs(ky + G_parameter_y); 
plot(ky, E, 'Linewidth', 2.5);
hold on 
plot(ky, sort(E, 'descend'), 'Linewidth', 2.5);

%******************
%% TE Mode 
%****************** 

E_TM = ((hev*c)/n).*sqrt(abs(ky.^2 + G_parameter_x^2)); 
plot(ky, E_TM, 'LineWidth', 2.5); 
xlabel('Ky'); ylabel('E(eV)'); title('Rectangular Lattice Dispersion relation', 'FontSize', 14); 
legend('TE(0,1)','TE (0,-1)', 'TM'); 
colorbar
