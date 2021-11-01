function [U_AC_V, F_AC_V, F_AC] = Gorkov(P_avg,v_avg)
%GORKOV Summary of this function goes here
%   Detailed explanation goes here

Parameters;

U_AC_V=(f_2/(2*rho_oil*v_0oil^2)*P_avg-f_2*(3/4)*rho_oil*v_avg); % Her udregnes gorkovs potential pr volumen af partikel

F_AC_V=-diff(U_AC_V); % Kraften pr volumen.

F_AC=F_AC_V*V_particle; % Kraften på en partikel til et givet givet sted i bølgen
end

