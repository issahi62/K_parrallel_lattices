%% 2D PBG:Square lattice
%% sqband.m clear all
%warning off tic;
epsa=8.9; 
epsb=1;
a=1.0; 
R=0.2*a; 
i=sqrt(-1); 
f=pi*R^2/a^2; 
NumberofCell=5; 
a1=a;
a2=a*i;
b1=2*pi/a/NumberofCell;
b2=2*pi/a/NumberofCell*i;
n=input('please input n:');
disp('Fourier transforming.....');
NumberofPW=(2*n+1)^2; 
mind=(-NumberofPW:NumberofPW)'+NumberofPW+1;
mind=mind(:,ones(1,size(mind,1)))-NumberofPW-1; 
GG=mind'*b1+mind*b2;
%clear mind;
eps21=2*f*( epsa-epsb)*besselj(1,abs(GG).*R)./(abs(GG).*R); 
eps21(NumberofPW+1,NumberofPW+1)= epsb+f*(epsa-epsb); 
%zz=[0,0]*[a1 a2].';
%use 1X1 supercell to verify the algorithm %5X5 super cell
zz=[
-2 -2;-2 -1;-2 0;-2 1;-2 2;
-1 -2;-1 -1;-1 0;-1 1;-1 2;
0 -2; 0 -1; 0 1; 0,2;
1 -2; 1 -1; 1 0; 1 1; 1 2;
2,-2; 2,-1; 2 0; 2,1; 2,2]*[a1 a2].';
eps20=zeros(length(eps21));

for x=1:length(zz)
eps20=eps20+exp(i*(real(GG).*real(zz(x))+imag(GG).*imag(zz(x)))).*eps21;
end
ff=pi*R^2*length(zz)/(NumberofCell^2*a^2);
eps20=eps20./NumberofCell^2; 
eps20(NumberofPW+1,NumberofPW+1)= epsb+ff*(epsa-epsb);

%clear GG;
count=1; r=ones(2*n*NumberofCell+1,2); 
for y=-n:n
    for x=-n:n
        G(count)=x*b1+y*b2;
        r(count,:)=[ x,y];
        count=count+1;
    end
end
disp('Building eps(G,G) matrix from the FFT matrix');

for x=1:NumberofPW
    for y=x:NumberofPW
        b=r(x,:)-r(y,:)+(2*n+1)^2+1; 
        eps2(x,y)=eps20(b(1),b(2));
        eps2(y,x)=eps2( x,y);
    end
end

k1=2*pi/a*0.5.*(0:0.1:1);
k2=2*pi/a*(0.5+(0.1:0.1:1).*0.5*i);
k3=2*pi/a*(0.5+0.5*i).*(0.9:-0.1:0);
k0=[k1 k2 k3]./NumberofCell;
disp('Now calculate the eigen values..'); 
eps2=inv(eps2);
if max(max(real(eps2))) > 10^7*max(max(imag(eps2)))
    disp('Your lattice is inversion symmetric'); 
    eps2=real(eps2);
else
disp('Imaginary part of FFT is not zero');
%stop;
%here we only demonstrate the inversion symmetric
%case %%however, nonsymmetric case is also supported
end
counter=1;
for ii=1:length(k0)
k=k0(ii);
%k=k0;
%M=(real(k+G.')*real(k+G)+imag(k+G.')*imag(k+G)).*eps2; %TE wave 
M=abs(k+G.')*abs(k+G).*eps2; %TM wave
E =sort(abs( eig(M)));
freq(:,counter)= sqrt(abs(E)).*a./2./pi;
fprintf('calculation of k=%f+%fi is finished\n',real(k),imag(k)); counter=counter+1;
end
tmpx=1:length(k0);
plot(tmpx,freq,'linewidth',2)
title('TM:Band structure of a 2D square PBG with a point defect (5X5)')
xlabel('wave vector')
ylabel('wa/2\pic')
grid on
axis([1 31 0.3 0.5])
tmp=ifft2(ifftshift(eps20));
%surf(fftshift(real(tmp))),shading interp,view(2)
%%to get the lattice, use tmp=ifft2(ifftshift(eps20));