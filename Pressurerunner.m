clc
clear
Parameters;
f = 10^6;

omega = 2*pi*f;

z=0.2;

z=[r_transducer:0.01:0.5];
p_z=zeros(1,length(z));
for n=1:length(z)
[p_z(n),C]=Pressure(z(n),omega);

end

plot(z,)




