n=H_33*C_0;
v_0 = sqrt(c33D/rho_P); % wave speed of compressional waves in the piezoelectric plate
k=omega/v_0; % Wave number 
Z0a = rho_P*v_0*S; % plane wave acoustic impedance of the piezoelectric plate.
Zba = 0; % <- vides ikke endnu, men det er: Corresponding acoustic impedance of the backing.
TeA = [1/n n/1i*omega*C_0
-1i*omega*C_0 0];
TaA = 1/(Zba-1i*Z0a*tan(kd/2))*[Zba+1i*Z0a*cot(kd) (Z0a)^2 +1i*Z0a*Zba*cot(kd)
1 Zab-2i*Z0a*tan(kd/2)];

TA = (TeA)*(TaA);