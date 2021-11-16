function [diffvect] = differentialboy(x,stepsize)
%DIFFERENTIALBOY Summary of this function goes here
%   Detailed explanation goes here

m=length(x);
for n=1:m
    if n==1 
      diffvect(n) = (x(n+1)-x(n)) / stepsize ;
    elseif n==m
       diffvect(n)=(x(n)-x(n-1)) / (stepsize);
    else
        diffvect(n) = (x(n+1)-x(n-1)) / (2*stepsize) ;
    end
end

