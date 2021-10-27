function [U_ac] = Gorkov(Wave1A, Wave2A)
%GORKOV Summary of this function goes here
%   Detailed explanation goes here

Parameters;
syms r z t
f_1 = 1 - kappa_p/kappa_l; % En konstant der skal bruges i gorkov. Burde nok være i parameters

f_2 = (2*(rho_p-rho_oil)) / (2*rho_p + rho_oil); % Samme som f_1

Wave1=Wave1A*sin(2*pi*z/lambda+omega*t); % Konstruerer den første trykbølge(venstre til højre). Burde måske være cos ifølge luke

Wave2=Wave2A*sin(2*pi*z/lambda-omega*t); % Konstruerer den anden trykbølge(højre til venstre). Amplituden er mindre end Wave1. Også måske cos

Wavesum=Wave1+Wave2; % Her lægges de to bølger sammen til et standing wave. Tvivlsomt tryk på 15 bar.

P_avg=f*int(Wavesum^2),t,0,(1/f); % Her udregnes det gennemsnitlige tryk over en periode <P>. press x to doubt matlab integration funktion

v=Wave1A/(rho_oil*v_0oil)*cos(2*pi*z/lambda+omega*t)+Wave2A/(rho_oil*v_0oil)*cos(2*pi*z/lambda-omega*t); % Her konstrueres hastigheden af mediet for den stående bølge. Burde være 90 grader ude af fase med trykket ifølge luke.

v_avg=f*int((v^2),t,0,(1/f)); % Her udregnes <v>. Doubtful

U_AC_V=(f_2/(2*rho_oil*v_0oil^2)*P_avg-f_2*(3/4)*rho_oil*v_avg); % Her udregnes gorkovs potential pr volumen af partikel

F_AC_V=-diff(U_AC_V); % Kraften pr volumen.

F_AC=F_AC_V*V_particle; % Kraften på en partikel til et givet givet sted i bølgen
end

