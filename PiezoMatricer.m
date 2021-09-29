clear 
clc
PiezoParameters;
% input Parameters:
V_in=300.


j = sqrt(-1);
v_0 = sqrt(c33D/rho_P); % wave speed of compressional waves in the piezoelectric plate
v_0b = sqrt(ElasticModolusBacking/rho_b);

f = 1.055*10^6; %frekvens

omega = 2*pi*f; %Vinkelhastighed
n=h_33*C_0; % Matrix parameter
k=omega/v_0; % Wave number 
Z0a = rho_P*v_0*S; % plane wave acoustic impedance of the piezoelectric plate.
Zba = rho_b*v_0b*S; % <- vides ikke endnu, men det er: Corresponding acoustic impedance of the backing Z_b^A Vi gætter foreløbigt på =rho_backing*v_0,backing*S Problemet er at der står den er afhængig af frekvens. Plus den er cone-shaped.
TeA = [1/n n/j*omega*C_0;-j*omega*C_0 0];
TaA = 1/(Zba-j*Z0a*tan(k*d/2))*([Zba+j*Z0a*cot(k*d) (Z0a)^2+j*Z0a*Zba*cot(k*d)
1 Zba-2*j*Z0a*tan(k*d/2)]);
TA = (TeA)*(TaA); % Totale transfer matrix
ZrAa =v_sOil*rho_oil*S; % Acoustic impedance of radiating medium. In this case the oil
SvIA = 1/(ZrAa*TA(2,1)+TA(2,2)); % v_t/I_in 
ZinAe = (ZrAa*TA(1,1)+TA(1,2))/(ZrAa*TA(2,1)+TA(2,2)); % V_in/I_in
S_VF = ZrAa*SvIA/ZinAe;   % F/V_in whoohoooo
F=V_in*S_VF % Force from input voltage.






