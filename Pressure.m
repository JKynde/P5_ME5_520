function [p_z,C,Wave] = Pressure(z,omega,t)
%PREASSURE Summary of this function goes here
%   Detailed explanation goes here
%Defraction coefficient
Parameters;
k=omega/v_0Oil;
C = 1-exp(j*k*(r_transducer)^2/(2*z));
[F,v_t] = Matricer(omega/(2*pi),300,1);
p_z =(rho_oil*v_0Oil*abs(v_t)*exp(j*k*z)); % Lige nu ikke ganget med C

R_reflect = (-rho_f*v_0f+rho_oil*v_0Oil) / (rho_oil*v_0Oil + rho_f*v_0f);

lambda = v_0Oil/10^6;
Wave = (2.2569*10^5/(1-R_reflect^2))*sin( (2*pi/lambda)*z + omega*t) + (2.2569*10^5/(R_reflect-R_reflect^3))*sin( (2*pi/lambda)*z - omega*t);

%Wave = 2*2.2569e+05*sin(2*pi*z/lambda)*cos(omega*t);





end

