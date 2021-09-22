PiezoParameters.m;
n=h_33*C_0; % Matrix parameter
v_0 = sqrt(c33D/rho_P); % wave speed of compressional waves in the piezoelectric plate
k=omega/v_0; % Wave number 
Z0a = rho_P*v_0*S; % plane wave acoustic impedance of the piezoelectric plate.
Zba = 0; % <- vides ikke endnu, men det er: Corresponding acoustic impedance of the backing Z_b^A(omega)
TeA = [1/n n/1i*omega*C_0
-1i*omega*C_0 0];
TaA = 1/(Zba-1i*Z0a*tan(k*d/2))*[Zba+1i*Z0a*cot(k*d) (Z0a)^2 +1i*Z0a*Zba*cot(k*d)
1 Zab-2i*Z0a*tan(k*d/2)];

TA = (TeA)*(TaA);

% Ligningesystemet: TA*[F;v]=[V*I] skal løses.