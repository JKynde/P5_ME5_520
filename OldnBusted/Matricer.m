function [F,v_t,ZinAe] = Matricer(f,V_in,param)
%Matricer Matricer bygger sittig matricerne og har Kraften og hastigheden
%ved transducerhovedet samt den elektriske impedans som output.
j = sqrt(-1);
omega = 2*pi*f; %Vinkelhastighed
k = omega/param.v_0; % Wave number piezo.
k_a = omega/param.v_0f;  % Wave number front layer
k_c = omega/param.v_0c;  % Wave number of cable
TeA = [1/param.n param.n/j*omega*param.C_0;-j*omega*param.C_0 0]; % Elektrisk matrice
TaA = 1/(param.Zba-j*param.Z0a*tan(k*param.d_p/2))*([param.Zba+j*param.Z0a*cot(k*param.d_p) (param.Z0a)^2+j*param.Z0a*param.Zba*cot(k*param.d_p)
1 param.Zba-2*j*param.Z0a*tan(k*param.d_p/2)]); % Acoustisk piezo matrice.
Tl = [cos(k_a*param.l_a) -j*param.Z0a_f*sin(k_a*param.l_a) % Acoustic front layer matrix
-j*sin(k_a*param.l_a)/param.Z0a_f cos(k_a*param.l_a)];% Front layer matrice
if param.Tlmode==1
    TA = ((TeA)*(TaA)*Tl); % Total transfer matrix
else 
    TA=(TeA)*(TaA);% Total transfer matrix
end
SvIA = 1/(param.ZrAa*TA(2,1)+TA(2,2)); % v_t/I_in 
ZinAe = (param.ZrAa*TA(1,1)+TA(1,2))/(param.ZrAa*TA(2,1)+TA(2,2)); % V_in/I_in
S_VF = param.ZrAa*SvIA/ZinAe;   % F/V_in whoohoooo
S_vV = SvIA / ZinAe ; % v_t/V_in Transducer speed over input voltage.
v_t = V_in * S_vV ; 
F = V_in*S_VF; % Force from input voltage.

end
%%
% % Cabel model 
% T_c=[cos(k_c*l_c) -j*Z0e*sin(k_c*l_c
% -j*sin(k_c*l_c)/Z0e cos(k_c*l_c)]
% %hvor
% Z0e=1/(2*pi)*sqrt(myh_c/epsilon_c)*ln(b_c/a_c);
% k_c = omega/v_0c;
