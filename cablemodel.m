function [Z] = cablemodel(f,V_in,param)
% https://en.wikipedia.org/wiki/Skin_effect
rho_m=1/param.sigma_c;
k=sqrt(-j*2*pi*f*param.mu_c/rho_m); %Complex wave number a donno man
R=param.b_c; % Ydre radius kabel

a=-pi;
b=pi;
t_steps=1000;
t_stepsize=(b-a)/(2*t_steps);
t=[a:t_stepsize:b-t_stepsize];
y_0=zeros(1,length(t));
y_1=zeros(1,length(t));
% Udregning af J_0
for n=1:length(t)
y_0(n)=exp(j*(k*R*sin(t(n))));
end
[j_0]=integralboy(y_0,t_stepsize);
J_0=1/(2*pi)*j_0;

% J_1(k*R)
for n=1:length(t)
y_1(n)=exp( j*( k*R*sin(t(n))-t(n) ) );
end
[j_1]=integralboy(y_1,t_stepsize);
J_1=1/(2*pi)*j_1;

ZprL=(k*rho_m)/(2*pi*R) * (J_0)/(J_1);
Z=ZprL*param.l_c*2;
end

