function [F,v_t,ZinAe] = Matricer2(f,V_in,Tlmode)
%Matricer Matricer bygger sittig matricerne og har Kraften og hastigheden
%ved transducerhovedet samt den elektriske impedans som output.
Parameters;
j = sqrt(-1);
omega = 2*pi*f; %Vinkelhastighed
k=omega/v_0; % Wave number piezo.
k_a = omega/v_0f;  % Wave number front layer
TeA = [1/n n/j*omega*C_0;-j*omega*C_0 0]; % Elektrisk matrice
TaA = 1/(Zba-j*Z0a*tan(k*d_p/2))*([Zba+j*Z0a*cot(k*d_p) (Z0a)^2+j*Z0a*Zba*cot(k*d_p)
1 Zba-2*j*Z0a*tan(k*d_p/2)]); % Acoustisk piezo matrice.
Tl = [cos(k_a*l_a) -j*Z0a_f*sin(k_a*l_a) % Acoustic front layer matrix
-j*sin(k_a*l_a)/Z0a_f cos(k_a*l_a)];% Front layer matrice
if Tlmode==1
TA = ((TeA)*(TaA)*Tl); % Total transfer matrix
else 
    TA=(TeA)*(TaA);% Total transfer matrix
end
SvIA = 1/(ZrAa*TA(2,1)+TA(2,2)); % v_t/I_in 
ZinAe = (ZrAa*TA(1,1)+TA(1,2))/(ZrAa*TA(2,1)+TA(2,2)); % V_in/I_in
S_VF = ZrAa*SvIA/ZinAe;   % F/V_in whoohoooo
S_vV = SvIA / ZinAe ; % v_t/V_in Transducer speed over input voltage.
v_t = V_in * S_vV ; 
F = V_in*S_VF; % Force from input voltage
end

