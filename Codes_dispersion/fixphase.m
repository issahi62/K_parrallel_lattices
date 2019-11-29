function z=fixphase(f)
s0=sum(sum(real(f).^2-imag(f).^2)); 
s1=sum(sum(2.*real(f).* imag(f))); 
theta=0.5*atan2(-s1,s0);
fprintf('Phase is:%f\n',theta); 
tmpr=real(f)* cos(theta)-imag(f)*sin(theta); 
tmpi=real(f)*sin(theta)-imag(f)*cos(theta);
z =tmpr+sqrt(-1)* tmpi;
end