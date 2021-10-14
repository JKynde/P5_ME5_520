function [U_ac] = Gorkov(p,v,r)
%GORKOV Summary of this function goes here
%   Detailed explanation goes here

Parameters;

syms r p
v = v_0Oil;
f_1 = 1 - kappa_p/kappa_l;

f_2 = (2*(rho_p-rho_oil)) / (2*rho_p + rho_oil);


U_ac = 4*pi/3 * r^3 * (f_1*1/(2*rho_oil*v_0Oil^2)*p^2 -f_2*3/4*rho_oil*v^2 )


end

