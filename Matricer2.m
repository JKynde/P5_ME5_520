function [F,v_t,ZinAe] = Matricer2(f,V_in,param)
%Matricer Matricer bygger sittig matricerne og har Kraften og hastigheden
%ved transducerhovedet samt den elektriske impedans som output.
Parameters;
j = sqrt(-1);
omega = 2*pi*f; %Vinkelhastighed
k=omega/param.v_0; % Wave number piezo.
k_a = omega/param.v_0f;  % Wave number front layer
N=1; % Remnent of old times. Must be 1, don't touch.
%% 4x4 matrice konstruktion og produkt.
theta = omega*param.d_p/param.v_0;
sigma = param.C_0*param.h_33^2/(omega*param.Z0a);
cosphi = (cos(theta)-sigma*sin(theta))/(1-sigma*sin(theta));
phi = acos(cosphi);
cosNphi = cos(N*phi);
sinNphi = sin(N*phi);
R = sqrt((sin(theta) - 2*sigma*(1-cos(theta)))/(sin(theta)));

T_disc = [cosNphi -j*param.Z0a*R*sinNphi -param.h_33*param.C_0*tan(1/2 * phi)*sinNphi 0
-j*(param.Z0a)^(-1)*R^(-1)*sinNphi cosNphi -j*param.h_33*C_0*param.Z0a^(-1)*R^(-1)*tan(1/2*phi)*(cosNphi-(-1)^N) 0
0 0 (-1)^N 0
-j*param.h_33*param.C_0*param.Z0a^(-1)*R^(-1)*tan(1/2*phi)*(cosNphi-(-1)^N) -param.h_33*param.C_0*tan(1/2 * phi)*sinNphi j*(N*(-1)^N)*(1+2*sigma*R^(-1)*tan(1/2*phi)+sigma*R^(-1)*tan(1/2*phi)*tan(1/2*phi)*sinNphi)*omega*param.C_0 (-1)^N];

if param.mode == 1
T = T_disc*T_disc;
elseif param.mode == 2

    theta_g = omega * param.l_m / param.v_0m;

    T_terminal = [cos(theta_g) j*param.Zma*sin(theta_g) 0 0
    j*sin(theta_g)/param.Zma cos(theta_g) 0 0
    0 0 1 0
    0 0 0 1];

T = T_terminal*T_disc*T_terminal*T_disc*T_terminal;
elseif param.mode == 3

    theta_g = omega * param.l_m / param.v_0m;
    theta_a = omega * param.l_aa / param.v_0a;

    T_terminal = [cos(theta_g) j*param.Zma*sin(theta_g) 0 0
    j*sin(theta_g)/param.Zma cos(theta_g) 0 0
    0 0 1 0
    0 0 0 1];

T_adhesive = [cos(theta_a) j*param.Zaa*sin(theta_a) 0 0
    j*sin(theta_a)/param.Zaa cos(theta_a) 0 0
    0 0 1 0
    0 0 0 1];

T = T_adhesive*T_terminal*T_adhesive*T_disc*T_adhesive*T_terminal*T_adhesive*T_disc*T_adhesive*T_terminal*T_adhesive;
elseif param.mode == 4
    theta_g = omega * param.l_m / param.v_0m;
    theta_a = omega * param.l_aa / param.v_0a;
    k_c = omega/param.c_c;
    delta_s = 1/sqrt(pi*f*param.mu_c*param.sigma_c);
    Z_c = (1+j)/(param.sigma_c*delta_s)*param.l_c/(2*pi*param.b_c);

    T_terminal = [cos(theta_g) j*param.Zma*sin(theta_g) 0 0
    j*sin(theta_g)/param.Zma cos(theta_g) 0 0
    0 0 1 0
    0 0 0 1];

T_adhesive = [cos(theta_a) j*param.Zaa*sin(theta_a) 0 0
    j*sin(theta_a)/param.Zaa cos(theta_a) 0 0
    0 0 1 0
    0 0 0 1];

T_cable = [1 0 0 0
    0 1 0 0
    0 0 cos(k_c*param.l_c) -j*Z_c*sin(k_c*param.l_c) 
    0 0 -j*sin(k_c*param.l_c)/Z_c cos(k_c*param.l_c)];

T = T_cable*T_adhesive*T_terminal*T_adhesive*T_disc*T_adhesive*T_terminal*T_adhesive*T_disc*T_adhesive*T_terminal*T_adhesive;
else
    T=T_disc;
end
if param.Tlmode==1
 Tl = [cos(k_a*param.l_a) -j*param.Z0a_f*sin(k_a*param.l_a) 0 0
-j*sin(k_a*param.l_a)/param.Z0a_f cos(k_a*param.l_a) 0 0
0 0 1 0
0 0 0 1];
 T=T*Tl;
end

%% Nedkogning og 2x2 produkteri


T_nedkogt = [T(3,1)-T(3,3)*((T(2,1)*param.Zba+T(1,1)))/(T(2,3)*Zba+T(1,3)) T(3,2)-T(3,3)*(T(2,2)*param.Zba+T(1,2))/(T(2,3)*Zba+T(1,3))
T(4,1)-T(4,3)*(T(2,1)*param.Zba+T(1,1))/(T(2,3)*Zba+T(1,3)) T(4,2)-T(4,3)*(T(2,2)*param.Zba+T(1,2))/(T(2,3)*param.Zba+T(1,3))];
 %if param.Tlmode==1
 %Tl = [cos(k_a*param.l_a) -j*param.Z0a_f*sin(k_a*param.l_a)
%-j*sin(k_a*param.l_a)/param.Z0a_f cos(k_a*param.l_a)];
 %TA = (T_nedkogt*Tl); % Total transfer matrix
 %else 
 TA=T_nedkogt; % Total transfer matrix
 %end
 %% Sensitivity functions and outputs
SvIA = 1/(param.ZrAa*TA(2,1)+TA(2,2)); % v_t/I_in 
ZinAe = (param.ZrAa*TA(1,1)+TA(1,2))/(param.ZrAa*TA(2,1)+TA(2,2)); % V_in/I_in
S_VF = param.ZrAa*SvIA/ZinAe;   % F/V_in whoohoooo
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
 %why no push
 %pls now
