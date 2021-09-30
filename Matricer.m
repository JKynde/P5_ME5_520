function [F,TA] = Matricer(f,V_in,Tlmode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Parameters;
j = sqrt(-1);
omega = 2*pi*f; %Vinkelhastighed
k=omega/v_0; % Wave number piezo.
k_a = omega/v_0f;  % Wave number front layer
TeA = [1/n n/j*omega*C_0;-j*omega*C_0 0]; % Elektrisk matrice
TaA = 1/(Zba-j*Z0a*tan(k*d_p/2))*([Zba+j*Z0a*cot(k*d_p) (Z0a)^2+j*Z0a*Zba*cot(k*d_p)
1 Zba-2*j*Z0a*tan(k*d_p/2)]); % Acoustic piezo matrix.
Tl = [cos(k_a*l_a) -j*Z0a_f*sin(k_a*l_a) % Acoustic front layer matrix
-j*sin(k_a*l_a)/Z0a_f cos(k_a*l_a)];% Front layer matrice
if Tlmode==1
TA = (TeA)*(TaA)*Tl; % Total transfer matrix
else 
    TA=(TeA)*(TaA);% Total transfer matrix
end
SvIA = 1/(ZrAa*TA(2,1)+TA(2,2)); % v_t/I_in 
ZinAe = (ZrAa*TA(1,1)+TA(1,2))/(ZrAa*TA(2,1)+TA(2,2)); % V_in/I_in
S_VF = ZrAa*SvIA/ZinAe;   % F/V_in whoohoooo
F = V_in*S_VF; % Force from input voltage.
end

