clear 
clc
PiezoParameters;
j = sqrt(-1);
f = 1*10^6; %frekvens
omega = 2*pi*f; %Vinkelhastighed
n=h_33*C_0; % Matrix parameter
v_0 = sqrt(c33D/rho_P); % wave speed of compressional waves in the piezoelectric plate
v_0b = sqrt(ElasticModolusBacking/rho_b);
k=omega/v_0; % Wave number 
Z0a = rho_P*v_0*S; % plane wave acoustic impedance of the piezoelectric plate.
Zba = rho_b*v_0b*S; % <- vides ikke endnu, men det er: Corresponding acoustic impedance of the backing Z_b^A Vi gætter foreløbigt på =rho_backing*v_0,backing*S Problemet er at der står den er afhængig af frekvens. Plus den er cone-shaped.
TeA = [1/n n/j*omega*C_0;-j*omega*C_0 0];
TaA = 1/(Zba-j*Z0a*tan(k*d/2))*([Zba+j*Z0a*cot(k*d) (Z0a)^2+j*Z0a*Zba*cot(k*d)
1 Zba-2*j*Z0a*tan(k*d/2)]);
TA = (TeA)*(TaA); % Totale transfer matrix

% Ligningesystemet: TA*[F;v]=[V*I] skal løses.

%%%%%%%%%%%%%%%%%%%%%%Nedenstående er fra impedans modellering artiklen.
epsilon33S=C_0/(S*d); 
K = sqrt(h_33*epsilon33S); %Static electromechanical coupling factor. Brugt til impedansen set fra tranducer terminalerne
X_c0 = -1*j/(omega*C_0); % En eller anden modstand. Bruges til Z_electrical
v_a = 1500 ; % Wave speed of load TBD. 
k_a = omega/v_a ; % fuck mig
l_a = 0.1; % Væskesøjle. Fuck mig.
Z_c = Z0a; % Characteristic impedance, antages at være ens med Z0a for nu.
Z_a = 1500*2000*S; % Characteristic impedance of the load. Sættes til nul for nu. Skal nok udregnes for olien.
Z_m1 = X_c0*(K^2*(tan(k*d)/(k*d/2))*((Z_c*tan(k*d/2))/(Z_c*tan(k*d)+Z_a*tan(k_a*l_a))));
Z_m2 = X_c0*(K^2*(tan(k*d)/(k*d/2))*((Z_a/2)*tan(k_a*l_a)/(Z_c*tan(k*d)+Z_a*tan(k_a*l_a))));
Z_e = X_c0+Z_m1+Z_m2; % equivalent electrical impedance seen from transducer terminals.


