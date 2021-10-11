function [outputArg1,outputArg2] = Gorkov(p,v)
%GORKOV Summary of this function goes here
%   Detailed explanation goes here

v_rms = (1/(sqrt(2)))*v;
p_rms = (1/(sqrt(2)))*p;

f_1 = 1 - kappa_p/kappa_l;

f_2 = (2*(rho_p-rho_l)) / (2*rho_p + rho_l);

U_ac = 4*pi/3 * r^3 * (f_1*1/(2*rho_l*c_l^2)*p_rms -f_2*3/4*rho_l*v_rms );


end

