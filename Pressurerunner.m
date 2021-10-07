clc
clear
Parameters;
f = 10^6;
omega = 2*pi*f;

t=[0:1*10^-7:(pi/omega)+2*10^-7];
z=[r_transducer*2:0.0001:0.02]; 
p_z=zeros(length(t),length(z));
C=zeros(1,length(z));

for i=1:length(t)
    
for n=1:length(z)
[p_z(i,n),C(n)]=Pressure(z(n),omega,t(i));
end
Progress=(i/length(t))*100;
Progress
end

surf(z,t,real(p_z))









