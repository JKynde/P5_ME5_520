function [F,TA] = Matricer(f,V_in)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Parameters;
j = sqrt(-1);
omega = 2*pi*f; %Vinkelhastighed
k=omega/v_0; % Wave number 
TeA = [1/n n/j*omega*C_0;-j*omega*C_0 0]; % Elektrisk matrice
TaA = 1/(Zba-j*Z0a*tan(k*d_p/2))*([Zba+j*Z0a*cot(k*d_p) (Z0a)^2+j*Z0a*Zba*cot(k*d_p)
1 Zba-2*j*Z0a*tan(k*d_p/2)]); % Acoustic piezo matrice.
TA = (TeA)*(TaA); % Totale transfer matrix
SvIA = 1/(ZrAa*TA(2,1)+TA(2,2)); % v_t/I_in 
ZinAe = (ZrAa*TA(1,1)+TA(1,2))/(ZrAa*TA(2,1)+TA(2,2)); % V_in/I_in
S_VF = ZrAa*SvIA/ZinAe;   % F/V_in whoohoooo
F = V_in*S_VF; % Force from input voltage.
end

