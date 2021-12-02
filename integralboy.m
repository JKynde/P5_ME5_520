function [Ettal] = integralboy(x,stepsize)
%INT Summary of this function goes here
%   Detailed explanation goes her
m=length(x);
Ettal=0;
for n=1:m
if n==1 | n==m
    Ettal=Ettal+x(n);
elseif mod(n,2)==0
    Ettal=Ettal+4*x(n);
elseif mod(n,2)==1
    Ettal=Ettal+2*x(n);
else end
end
Ettal=Ettal*stepsize/3;
end

