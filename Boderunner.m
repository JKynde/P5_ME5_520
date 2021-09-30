% Det her script kører Matricer.m funktionen med det formål at lave et
% bodeplot
clc
clear
j = sqrt(-1);
Parameters; % Først defineres alle parametre fra konstante Parameters.m
V_in = 300*exp(j*deg2rad(0)); % Indgangs spændingsvisor.

f = logspace(5,6.5,500); % Lav en vektor logaritmisk fordelte indgange med n punkter.
F_out = zeros(1,length(f)); % Lav en vektor med nuller lige så lang som vektor f.

for n=1:length(f) % For loop som fylder F_out med resultaterne af modellen
    
    F_out(n)=Matricer(f(n),V_in);
   
end
F_DB = 20*log10(abs(F_out));
figure;
semilogx(f,F_DB)
figure;
semilogx(f,rad2deg(angle(F_out)))

