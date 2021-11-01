%Dette script bruges til at vurdere om vores numeriske integration virker
clc
clear all
b=10;
a=-10;
t_steps=5000;
t_stepsize=(b-a)/(2*t_steps);
t=[a:t_stepsize:b-t_stepsize];
for n=1:length(t)
y(n)=4*t(n)^3+exp(t(n))+cos(t(n));
end
[ettal]=integralboy(y,t_stepsize)

envector = differentialboy(y,t_stepsize);

plot(t,envector)