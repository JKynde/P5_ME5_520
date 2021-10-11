clc
clear
Parameters;
f = 10^6;
omega = 2*pi*f;
t=[0:1*10^-8*2:2*10^-6];
z=[r_transducer*2:0.0001:2*r_transducer+0.002]; 
p_z=zeros(length(t),length(z));
Wave = zeros(length(t),length(z));
C=zeros(1,length(z));

for i=1:length(t)
    
for n=1:length(z)
[p_z(i,n),C(n),Wave(i,n)]=Pressure(z(n),omega,t(i));
end
Progress=(i/length(t))*100;
Progress
end



for i=1:length(t)
    hold on
    plot(z,Wave(i,:))
    pause(0.1)
    
end
    
    
    
    









