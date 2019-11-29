clear all
clc
%function omega = dispersion(phi)
syms L N omega k kx ky a b K G M
%Test Input
L=10;
N=10;
a=4.33;
b=20;
K=1;
G=2;
M=1;
j=1;
%Set loop to cover entire lattice
e=-N/2;
r=N/2;
for i = e:.1:r
 %Wave vector
 kx(j)= (2*pi()/a)*(i/L);
 %Dispersion Curve Equation
 omega1(j) = sqrt(((K+G)/M)+(1/M)*sqrt(K^2+G^2+2*K*G*cos(kx(j)*a)));
 omega2(j) = sqrt(((K+G)/M)-(1/M)*sqrt(K^2+G^2+2*K*G*cos(kx(j)*a)));
 j=j+1;
end
kx=double(kx);
omega1=double(omega1);
omega2=double(omega2);
plot(kx, omega1, kx, omega2)
legend('Optical Branch','Acoustic Branch','Location','NorthEast')
xlabel('k')
ylabel('w(k)')
axis([-pi()/a pi()/a 0 3.5]) 