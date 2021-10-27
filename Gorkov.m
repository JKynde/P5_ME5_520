function [U_ac] = Gorkov(Wave1A, Wave2A)
%GORKOV Summary of this function goes here
%   Detailed explanation goes here

Parameters;
V=(4/3)*pi*(10^-4)^3
syms r z t
f_1 = 1 - kappa_p/kappa_l;

f_2 = (2*(rho_p-rho_oil)) / (2*rho_p + rho_oil);

Wave1=Wave1A*sin(2*pi*z/lambda+omega*t)

Wave2=Wave2A*sin(2*pi*z/lambda-omega*t)

Wavesum=Wave1+Wave2

P_avg=f*int(Wavesum^2),t,0,(1/f)

v=Wave1A/(rho_oil*v_0oil)*cos(2*pi*z/lambda+omega*t)+Wave2A/(rho_oil*v_0oil)*cos(2*pi*z/lambda-omega*t)

v_avg=f*int((v^2),t,0,(1/f))

U_AC_V=(f_2/(2*rho_oil*v_0oil^2)*P_avg-f_2*(3/4)*rho_oil*v_avg)

F_AC_V=-diff(U_AC_V)

F_AC=F_AC_V*V
end

