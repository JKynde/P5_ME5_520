clear 
clc
PiezoParameters;
omega = 1; %frekvens?
n=h_33*C_0; % Matrix parameter
v_0 = sqrt(c33D/rho_P); % wave speed of compressional waves in the piezoelectric plate
k=omega/v_0; % Wave number 
Z0a = rho_P*v_0*S; % plane wave acoustic impedance of the piezoelectric plate.
Zba = 1; % <- vides ikke endnu, men det er: Corresponding acoustic impedance of the backing Z_b^A(omega)
TeA = [1/n n/1i*omega*C_0;-1i*omega*C_0 0];
TaA = 1/(Zba-1i*Z0a*tan(k*d/2))*([Zba+1i*Z0a*cot(k*d) (Z0a)^2+1i*Z0a*Zba*cot(k*d)
1 Zba-2i*Z0a*tan(k*d/2)]);
TA = (TeA)*(TaA); % Totale transfer matrix

% Ligningesystemet: TA*[F;v]=[V*I] skal løses.

%%%%%%%%%%%%%%%%%%%%%%Nedenstående er fra impedans modellering artiklen.
epsilon33S=C_0/(areal*l); 
K = sqrt(h_33*epsilon33S); %Static electromechanical coupling factor. Brugt til impedansen set fra tranducer terminalerne
X_c0 = -1*j/(omega*C_0);
Z_a=rho_a*v_a*A; 
Z_m1 = X_c0*(K^2*(tan(k*l)/(k*l/2))*((Z_c*tan(k*l/2))/(Z_c*tan(k*l)+Z_a*tan(k_a*l_a))));
Z_m2 = X_c0*(K^2*(tan(k*l)/(k*l/2))*((Z_a/2)*tan(k_a*l_a)/(Z_c*tan(k*l)+Z_a*tan(k_a*l_a))));
Zelectrical = X_C0+Z_m1+Z_m2; % equivalent electrical impedance seen from transducer terminals.

