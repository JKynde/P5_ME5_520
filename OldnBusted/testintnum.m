%Dette script bruges til at vurdere om vores numeriske integration virker
clc
clear all
a=0;
b=2*pi;
t_steps=5000;
t_stepsize=(b-a)/(2*t_steps);
t=[a:t_stepsize:b-t_stepsize];
for n=1:length(t)
y(n)=sin(t(n));
end
[ettal]=integralboy(y,t_stepsize)

envector = differentialboy(y,t_stepsize);

plot(t,envector)