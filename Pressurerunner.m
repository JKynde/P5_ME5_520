clc
clear
Parameters;
f = 10^6;

omega = 2*pi*f;

z=[r_transducer*2:0.00001:0.05];
p_z=zeros(1,length(z));
C=zeros(1,length(z));
for n=1:length(z)
[p_z(n),C(n)]=Pressure(z(n),omega);
end

plot(z,p_z)





