function [p_z,C,WavesumP, F_z,Wave1P,Wave2P,Wavesumv] = Pressure(z,omega,t,F,v_t)
%PREASSURE Summary of this function goes here
%   Detailed explanation goes here
%Defraction coefficient
Parameters;
k=omega/v_0Oil; %Wave number for oil

C = 1-exp(j*k*(r_transducer)^2/(2*z)); % Diffraction coefficient anvendes ikke i øjeblikket. Fra p5biblen.

p_z =(rho_oil*v_0Oil*abs(v_t)*exp(j*k*z)); % Lige nu ikke ganget med C, men ellers er det paraxial approximation fra P5biblen.

P_tmax = abs(F)/Area;

Wave1A =(P_tmax/(1-R_reflect^2)); % Her udregnes amplituderne for de to bølger der udgør den stående bølge
Wave2A = (P_tmax/(R_reflect*(1-R_reflect^2))); % Deres amplituder regnes som en uende geometrisk række, som bliver ved med at blive reflekteret.

Wave1P = Wave1A*sin( (2*pi/lambda)*z + omega*t); % Bølgerne laves.
Wave2P = Wave2A*sin( (2*pi/lambda)*z - omega*t); 

Wave1v = Wave1A/(rho_oil*v_0Oil) * sin ( (2*pi/lambda)*z + omega*t - pi/2);
Wave2v = Wave2A/(rho_oil*v_0Oil) * sin ( (2*pi/lambda)*z - omega*t - pi/2);

WavesumP = Wave1P + Wave2P; %Standing wave equation. Summen af alle venstregående og højregående bølger.
Wavesumv = Wave1v + Wave2v;
%Wave = 2*2.2569e+05*sin(2*pi*z/lambda)*cos(omega*t);

% Lukes approximeret metode hvor amplituden af den stående bølge er
% 2*P_tmax. Desuden er den stående bølge også perfekt.
I_sound =(P_tmax)^2/(2*rho_oil*v_0Oil);
P_effective=I_sound*Area; % https://courses.lumenlearning.com/physics/chapter/17-3-sound-intensity-and-sound-level/ ikke lukas men link

F_z=((pi*(omega/(2*pi))*P_effective * (rho_oil + 11*rho_p)) / (2*Area*v_0Oil^2*(rho_oil + 2*rho_p)))*sin(2*pi*(omega/(2*pi))*z/v_0Oil);



end

