function Dispersion2D(epsilon, sigma, r, m)
%Phonon Dispersion Curves for Two-Dimensional Monatomic Lattice
%function Dispersion2D(epsilon, sigma, r, m)
%Calcualte Interaction Energy
phi= 8*epsilon*(((78*sigma^12)/(r^14))-((21*sigma^6)/(r^8)))...
+8*epsilon*(((78*sigma^12)/((2*r)^14))-((21*sigma^6)/((2*r)^8)));
%Next nearest neighbors addition
%Xenon Test Inputs
%epsilon=0.02;
%sigma=3.98;
%r=4.33;
%m=131.29
%sym L N omega k kx ky a b K G M
%Input for Monatomic Homogeneous Lattice:
a=r;
b=r;
K=phi;
G=phi;
M=m*1.66053886e-27; %convert amu to kg
%Lattice Size:
L=10; N=10;
%Set loop to cover entire lattice
e=N/2;
j=1;
for i = -e:.1:e
 kx(j)= (2*pi()/a)*(i/L);
 ky(j)= (2*pi()/b)*(i/N);
 %k(i)=[kx(i), ky(i)];
 %Dispersion Curve Equation
 omega(j) = 2*sqrt(K/M)*abs(sin((1/2)*kx(j)*a))...
 +2*sqrt(G/M)*abs(sin((1/2)*ky(j)*b));
 j=j+1;
end
kx=double(kx);
ky=double(ky);
omega=double(omega);
plot(kx, omega)
xlabel('k')
ylabel('w(k)')
axis([-pi()/a pi()/a 0 1.1*max(omega)])
legend('Acoustic Branch','Location','North')
figure
plot3(kx, ky, omega)
xlabel('k'); ylabel('k'); zlabel('w(k)') 