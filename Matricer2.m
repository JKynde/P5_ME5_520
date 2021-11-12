function [F,v_t,ZinAe] = Matricer2(f,V_in,Tlmode,N)
%Matricer Matricer bygger sittig matricerne og har Kraften og hastigheden
%ved transducerhovedet samt den elektriske impedans som output.
Parameters;
j = sqrt(-1);
omega = 2*pi*f; %Vinkelhastighed
k=omega/v_0; % Wave number piezo.
k_a = omega/v_0f;  % Wave number front layer
%% Sindsyg matrice fra sindsyg kilde:
theta = omega*d_p/v_0;
sigma = C_0*h_33^2/(omega*Z0a);
cosphi = (cos(theta)-sigma*sin(theta))/(1-sigma*sin(theta));
phi = acos(cosphi);
cosNphi = cos(N*phi);
sinNphi = sin(N*phi);

R = sqrt((sin(theta) - 2*sigma*(1-cos(theta)))/(sin(theta)));

T = [cosNphi -j*Z0a*R*sinNphi -h_33*C_0*tan(1/2 * phi)*sinNphi 0
-j*(Z0a)^(-1)*R^(-1)*sinNphi cosNphi -j*h_33*C_0*Z0a^-1*R^-1*tan(1/2*phi)*(cosNphi-(-1)^N) 0
0 0 (-1)^N 0
-j*h_33*C_0*Z0a^-1*R^-1*tan(1/2*phi)*(cosNphi-(-1)^N) -h_33*C_0*tan(1/2 * phi)*sinNphi j*(N*(-1)^N)*(1+2*sigma*R^-1*tan(1/2*phi)+sigma*R^(-1)*tan(1/2*phi)*tan(1/2*phi)*sinNphi) (-1)^N];





%% Hello




T_nedkogt = [T(3,1)-T(3,3)*((T(2,1)*Zba+T(1,1)))/(T(2,3)*Zba+T(1,3)) T(3,2)-T(3,3)*(T(2,2)*Zba+T(1,2))/(T(2,3)*Zba+T(3,1))
T(4,1)-(T(4,3)*(T(2,1)*Zba+T(1,1)))/(T(2,3)*Zba+T(3,1)) T(4,2)-T(4,3)*(T(2,2)*Zba+T(2,1))/(T(2,3)*Zba+T(3,1))];
Tl = [cos(k_a*l_a) -j*Z0a_f*sin(k_a*l_a)
-j*sin(k_a*l_a)/Z0a_f cos(k_a*l_a)];
 if Tlmode==1
 TA = (T_nedkogt*Tl); % Total transfer matrix
 else 
     TA=T_nedkogt;% Total transfer matrix
end
SvIA = 1/(ZrAa*TA(2,1)+TA(2,2)); % v_t/I_in 
ZinAe = (ZrAa*TA(1,1)+TA(1,2))/(ZrAa*TA(2,1)+TA(2,2)); % V_in/I_in
S_VF = ZrAa*SvIA/ZinAe;   % F/V_in whoohoooo
S_vV = SvIA / ZinAe ; % v_t/V_in Transducer speed over input voltage.
v_t = V_in * S_vV ; 
F = V_in*S_VF; % Force from input voltage
end
%k_m = omega/v_0m; % Wave number mellem layer
% f_piezo = v_0/(2*d_p); % f_mn fra sittig til piezolagene
% gamma_piezo = k_a*d_p;
% phi_piezo = k_33*sqrt(2*pi*f_piezo*C_0*Z0a/pi);
% s_sittig = k_33^2*sin(gamma_piezo)/gamma_piezo;
% c_sittig = (k_33^2*(1-cos(gamma_piezo)))/gamma_piezo;
% gamma_mellem = k_m*l_m; 

% T_24 = [(cos(gamma_piezo)-s_sittig)/(1-s_sittig) (j*Z0a*(sin(gamma_piezo)-2*c_sittig))/(1-s_sittig) -((cos(gamma_piezo)-1)*phi_piezo)/(1-s_sittig) 0
% j*sin(gamma_piezo)/(Z0a*(1-s_sittig)) (cos(gamma_piezo)-s_sittig)/(1-s_sittig) (-j*phi_piezo*sin(gamma_piezo))/(1-s_sittig) 0
% 0 0 1 0
% (-j*sin(gamma_piezo)*phi_piezo)/(Z0a*(1-s_sittig)) -(cos(gamma_piezo)-1)*phi_piezo/(1-s_sittig) j*omega*C_0/(1-s_sittig) 1];
% 
% % Matricer for elektriske mellemlag
% T_135 = [cos(gamma_mellem) j*Zma*sin(gamma_mellem) 0 0
% j*sin(gamma_mellem)/(Zma) cos(gamma_mellem) 0 0
% 0 0 1 0
% 0 0 0 1];
% 
% T = T_135*T_24*T_135*T_24*T_135;
% Vi skal have bygget 6 matricer i alt. Mellemlag elektrisk - Piezo -
% Mellemlag elektrisk - Piezo - Mellem ikke tilsluttet terminal og frontlag til sidst de 5 første koges ned til 2x2 og ganges
% med t_l som er front lag matricen som vi kender den.

%Matricer for piezodisc 1 og 2, hvilket er matrice T_2 og T_4 ergo T_24 .

