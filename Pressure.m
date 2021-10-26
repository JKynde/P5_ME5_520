function [p_z,C,Wavesum,F_z,Wave1,Wave2] = Pressure(z,omega,t)
%PREASSURE Summary of this function goes here
%   Detailed explanation goes here
%Defraction coefficient
Parameters;
k=omega/v_0Oil; %Wave number for oil

C = 1-exp(j*k*(r_transducer)^2/(2*z)); % Diffraction coefficient

[F,v_t] = Matricer(omega/(2*pi),300,1); % Force and velocity from Sittig Model
p_z =(rho_oil*v_0Oil*abs(v_t)*exp(j*k*z)); % Lige nu ikke ganget med C

R_reflect = (-rho_f*v_0f+rho_oil*v_0Oil) / (rho_oil*v_0Oil + rho_f*v_0f); %Reflection coefficient

Wave1A =(2.24701*10^5/(1-R_reflect^2)); % Her udregnes amplituderne for de to bølger der udgør den stående bølge
Wave2A = (2.24701*10^5/(R_reflect*(1-R_reflect^2))); % Deres amplituder regnes som en uende geometrisk række, som bliver ved med at blive reflekteret.

Wave1 = Wave1A*sin( (2*pi/lambda)*z + omega*t); % Bølgerne laves.
Wave2 = Wave2A*sin( (2*pi/lambda)*z - omega*t);

Wavesum = Wave1 + Wave2; %Standing wave equation. Summen af alle venstregående og højregående bølger.

%Wave = 2*2.2569e+05*sin(2*pi*z/lambda)*cos(omega*t);

I_sound =(2.24701*10^5)^2/(2*rho_oil*v_0Oil);
P_effective=I_sound*Area; % https://courses.lumenlearning.com/physics/chapter/17-3-sound-intensity-and-sound-level/

F_z=((pi*(omega/(2*pi))*P_effective * (rho_oil + 11*rho_p)) / (2*Area*v_0Oil^2*(rho_oil + 2*rho_p)))*sin(2*pi*(omega/(2*pi))*z/v_0Oil);



end

