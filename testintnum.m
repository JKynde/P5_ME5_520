%Dette script bruges til at vurdere om vores numeriske integration virker
clc
clear all
b=10;
a=0;
t_steps=5000;
t_stepsize=(b-a)/(2*t_steps);
t=[a:t_stepsize:b-t_stepsize];
for n=1:length(t)
y(n)=(3*sin(t(n))+4*cos(2*t(n)))^2;
end
[ettal]=integralboy(y,t_stepsize)