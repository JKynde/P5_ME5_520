function [p_z,C] = Pressure(z,omega,t)
%PREASSURE Summary of this function goes here
%   Detailed explanation goes here
%Defraction coefficient
Parameters;
k=omega/v_0Oil;
C = 1-exp(j*k*(r_transducer)^2/(2*z));
[F,v_t] = Matricer(omega/(2*pi),300,1);
p_z = C*(rho_oil*v_0Oil*abs(v_t)*exp(j*k*z));

%*exp(i*omega*t)


end

