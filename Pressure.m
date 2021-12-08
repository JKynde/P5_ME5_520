function [p_z,C,WavesumP, F_z,Wave1P,Wave2P,Wavesumv] = Pressure(z,omega,t,F,v_t,param)
%PRESSURE Pressure bygger trykfeltet og hastighedsfeltet ud fra en z
%vektor, vinkelhastigheden, en tidsvektor samt det maksimale tryk og
%hastighed ved transducerhovedet.
%Defraction coefficient
%Parameters;
k=omega/param.v_0Oil; %Wave number for oil

C = 1-exp(j*k*(param.r_transducer)^2/(2*z)); % Diffraction coefficient anvendes ikke i øjeblikket. Fra p5biblen.

p_z =(param.rho_oil*param.v_0Oil*abs(v_t)*exp(j*k*z)); % Lige nu ikke ganget med C, men ellers er det paraxial approximation fra P5biblen.

P_tmax = abs(F)/param.Area;

Vaac = (2*param.nu_oil*omega^2)/(3*(param.v_0Oil^3));% Viscous acoustic attenuation coefficient.
Vaac_one_trip=exp(-Vaac*param.d_heads); % Acoustic attenuation af a result of one trip between the heads.
Decay_constant = param.R_reflect*Vaac_one_trip; % Final decay constant pr trip for geometric series.

Wave1A = P_tmax/(1-Decay_constant^2); % Deres amplituder regnes som en uende geometrisk række, som bliver ved med at blive reflekteret.
Wave2A =Decay_constant*P_tmax/(1-Decay_constant^2); % Her udregnes amplituderne for de to bølger der udgør den stående bølge


Wave1P = Wave1A*exp(j*((k)*z+omega*t)); % Bølgerne laves som fasere
Wave2P = Wave2A*exp(j*((-k)*z+omega*t));

Wave1v = 1/(param.rho_oil*param.v_0Oil) * Wave1A*exp(j*((k)*z+omega*t));
Wave2v = -1/(param.rho_oil*param.v_0Oil) * Wave2A*exp(j*((-k)*z+omega*t));

WavesumP = real(Wave1P + Wave2P); %Standing wave equation. Summen af alle venstregående og højregående bølger.
Wavesumv = real(Wave1v + Wave2v);

% Lukes approximeret metode hvor amplituden af den stående bølge er
% 2*P_tmax. Desuden er den stående bølge også perfekt.
I_sound =(P_tmax)^2/(2*param.rho_oil*param.v_0Oil);
P_effective=I_sound*param.Area; % https://courses.lumenlearning.com/physics/chapter/17-3-sound-intensity-and-sound-level/ ikke lukas men link

F_z=((pi*(omega/(2*pi))*P_effective * (param.rho_oil + 11*param.rho_p)) / (2*param.Area*param.v_0Oil^2*(param.rho_oil + 2*param.rho_p)))*sin(2*pi*(omega/(2*pi))*z/param.v_0Oil);



end

