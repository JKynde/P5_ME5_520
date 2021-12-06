function [F,v_t,ZinAe,S_VF] = Matricer2(f,V_in,param)
%Matricer Matricer bygger sittig matricerne og har Kraften og hastigheden
%ved transducerhovedet samt den elektriske impedans som output.
j = sqrt(-1);
omega = 2*pi*f; %Vinkelhastighed
k=omega/param.v_0; % Wave number piezo.
k_a = omega/param.v_0f;  % Wave number front layer
N=1; % Remnent of old times. Must be 1, don't touch.
%% 4x4 matrice konstruktion og produkt.
theta = omega*param.d_p/param.v_0; %random parameter fra bloomfield
sigma = param.C_0*param.h_33^2/(omega*param.Z0a); %same
cosphi = (cos(theta)-sigma*sin(theta))/(1-sigma*sin(theta)); %cos(phi) bygges
phi = acos(cosphi); %phi defineres
cosNphi = cos(N*phi); %cos(Nphi) defineres
sinNphi = sin(N*phi); %sin(Nphi) defineres
R = sqrt((sin(theta) - 2*sigma*(1-cos(theta)))/(sin(theta))); %R defineres

% T_disc defineres ud fra bloomfield. Matrice er egenligt for flere ens
% disce.
T_disc = [cosNphi -j*param.Z0a*R*sinNphi -param.h_33*param.C_0*tan(1/2 * phi)*sinNphi 0
-j*(param.Z0a)^(-1)*R^(-1)*sinNphi cosNphi -j*param.h_33*param.C_0*param.Z0a^(-1)*R^(-1)*tan(1/2*phi)*(cosNphi-(-1)^N) 0
0 0 (-1)^N 0
-j*param.h_33*param.C_0*param.Z0a^(-1)*R^(-1)*tan(1/2*phi)*(cosNphi-(-1)^N) -param.h_33*param.C_0*tan(1/2 * phi)*sinNphi j*(N*(-1)^N)*(1+2*sigma*R^(-1)*tan(1/2*phi)+sigma*R^(-1)*tan(1/2*phi)*tan(1/2*phi)*sinNphi)*omega*param.C_0 (-1)^N];

if param.mode == 1 %Helt simpelt case hvor der kun er disce
T = T_disc*T_disc;
elseif param.mode == 2 %case med elektriske terminaler

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

T_adhesive = [cos(theta_a) -j*param.Zaa*sin(theta_a) 0 0
    -j*sin(theta_a)/param.Zaa cos(theta_a) 0 0
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

T_adhesive = [cos(theta_a) -j*param.Zaa*sin(theta_a) 0 0
    -j*sin(theta_a)/param.Zaa cos(theta_a) 0 0
    0 0 1 0
    0 0 0 1];
 
T_cable = [1 0 0 0
    0 1 0 0
    0 0 cos(k_c*param.l_c) -j*Z_c*sin(k_c*param.l_c) 
    0 0 -j*sin(k_c*param.l_c)/Z_c cos(k_c*param.l_c)];

% T_cable = [cos(k_c*param.l_c) -j*Z_c*sin(k_c*param.l_c) 0 0
% -j*sin(k_c*param.l_c)/Z_c cos(k_c*param.l_c) 0 0
% 0 0 1 0
% 0 0 0 1];

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


T_nedkogt = [T(3,1)-T(3,3)*((T(2,1)*param.Zba+T(1,1)))/(T(2,3)*param.Zba+T(1,3)) T(3,2)-T(3,3)*(T(2,2)*param.Zba+T(1,2))/(T(2,3)*param.Zba+T(1,3))
T(4,1)-T(4,3)*(T(2,1)*param.Zba+T(1,1))/(T(2,3)*param.Zba+T(1,3)) T(4,2)-T(4,3)*(T(2,2)*param.Zba+T(1,2))/(T(2,3)*param.Zba+T(1,3))];
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