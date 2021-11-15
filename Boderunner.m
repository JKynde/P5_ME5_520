% Det her script kører Matricer.m funktionen med det formål at lave et
% bodeplot
clc
clear
j = sqrt(-1);
plotmode = 1; % 1 = F   2 = Z    3 = v     4 = F2
Parameters; % Først defineres alle parametre fra konstante Parameters.m
V_in = 150*exp(j*deg2rad(0)); % Indgangs spændingsvisor.
Tlmode = 1; % Tlmode til Matricer.m functionen. Der er om front layer matricen er ganget på. 1 for ja else ikke.

f = logspace(5,6.5,5000); % Lav en vektor logaritmisk fordelte indgange med n punkter.
F_out = zeros(1,length(f)); % Lav en vektor med nuller lige så lang som vektor f.
Z_out = zeros(1,length(f)); 
v_out = zeros(1,length(f));
F_out2 = zeros(1,length(f));
for n=1:length(f) % For loop som fylder F_out med resultaterne af modellen
    %F_out(n)=Matricer(f(n),V_in,Tlmode);
    [F_out(n),v_out(n),Z_out(n)]=Matricer(f(n),V_in,Tlmode); %for loopet løber over forskellige værdier af f og outputter tre vektorer som er F, Z, V
    [F_out2(n)] = Matricer2(f(n),V_in,Tlmode,1); % det samme som ovenover men for sittig2
end

F_DB2 = 20*log10(abs(F_out2)); 
F_DB = 20*log10(abs(F_out)); % Her omregnes værdierne til db værdier
Z_DB = 20*log10(abs(Z_out));
v_DB = 20*log10(abs(v_out));


switch plotmode
    case 1 %F
        figure;

t=tiledlayout(2,1); % Her plottes outputtet af matricer som bodeplot
xlabel(t,'log spaced f values')

nexttile
semilogx(f,F_DB) % skriv enten F_DB, V_DB eller Z_DB
grid on
title('Magnitude plot of F(f)')
ylabel('Magnitude of F(dB)')

nexttile
semilogx(f,rad2deg(angle(F_out))) % Samme her
grid on
title('Angle plot of F(f)')
ylabel('Angle of F(^\circ)')
    case 2 %Z
        figure;

t=tiledlayout(2,1); % Her plottes outputtet af matricer som bodeplot
xlabel(t,'log spaced f values')

nexttile
semilogx(f,Z_DB) % skriv enten F_DB, V_DB eller Z_DB
grid on
title('Magnitude plot of F(f)')
ylabel('Magnitude of F(dB)')

nexttile
semilogx(f,rad2deg(angle(Z_out))) % Samme her
grid on
title('Angle plot of F(f)')
ylabel('Angle of F(^\circ)')
    case 3 %v
        figure;

t=tiledlayout(2,1); % Her plottes outputtet af matricer som bodeplot
xlabel(t,'log spaced f values')

nexttile
semilogx(f,v_DB) % skriv enten F_DB, V_DB eller Z_DB
grid on
title('Magnitude plot of F(f)')
ylabel('Magnitude of F(dB)')

nexttile
semilogx(f,rad2deg(angle(v_out))) % Samme her .
grid on
title('Angle plot of F(f)')
ylabel('Angle of F(^\circ)')

case 4 %F2
        figure;

t=tiledlayout(2,1); % Her plottes outputtet af matricer som bodeplot
xlabel(t,'log spaced f values')

nexttile
semilogx(f,F_DB2) % skriv enten F_DB, V_DB eller Z_DB
grid on
title('Magnitude plot of F2(f)')
ylabel('Magnitude of F2(dB)')

nexttile
semilogx(f,rad2deg(angle(F_out2))) % Samme her .
grid on
title('Angle plot of F2(f)')
ylabel('Angle of F2(^\circ)')
        
    otherwise
end


