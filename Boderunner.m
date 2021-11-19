% Det her script kører Matricer.m funktionen med det formål at lave et
% bodeplot
clc
clear all
close all
Parameters; % Først defineres alle parametre fra konstante Parameters.m
j = sqrt(-1);
plotmodes = 2; % 1 = F   2 = Z    3 = v     4 = F Matricer2   5 = Z Matricer2
V_in = 150*exp(j*deg2rad(0)); % Indgangs spændingsvisor.
mode = 3; % Running modes for Matricer2 1 = 2 diske 2 = 2 diske + terminaler 3 = 2 disk + terminaler + lim else 1 disk ligesom matricer
Tlmode = 1; % Tlmode til Matricer.m functionen. Der er om front layer matricen er ganget på. 1 for ja else ikke.
f = logspace(5,6.5,5000); % Lav en vektor logaritmisk fordelte indgange med n punkter.
F_out = zeros(1,length(f)); % Lav en vektor med nuller lige så lang som vektor f.
Z_out = zeros(1,length(f)); 
v_out = zeros(1,length(f));
for n=1:length(f) % For loop som fylder F_out med resultaterne af modellen
    if plotmodes == 1 || plotmodes == 2 || plotmodes == 3
        [F_out(n),v_out(n),Z_out(n)]=Matricer(f(n),V_in,Tlmode); %for loopet løber over forskellige værdier af f og outputter tre vektorer som er F, Z, V
    elseif plotmodes == 4 || plotmodes == 5  
        [F_out(n),Z_out(n)] = Matricer2(f(n),V_in,Tlmode,mode); % det samme som ovenover men for sittig2   
    else
    end
 end

%F_DB = 20*log10(abs(F_out2)); 
%Z_in_DB = 20*log10(abs(Z_in));
F_DB = 20*log10(abs(F_out)); % Her omregnes værdierne til db værdier
Z_DB = 20*log10(abs(Z_out));
v_DB = 20*log10(abs(v_out));


switch plotmodes
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

case 4 %F med sittig/bloomfield
        figure;

t=tiledlayout(2,1); % Her plottes outputtet af matricer som bodeplot
xlabel(t,'log spaced f values')

nexttile
semilogx(f,F_DB) % skriv enten F_DB, V_DB eller Z_DB
grid on
title('Magnitude plot of F(f)')
ylabel('Magnitude of F(dB)')

nexttile
semilogx(f,rad2deg(angle(F_out))) % Samme her .
grid on
title('Angle plot of F(f)')
ylabel('Angle of F(^\circ)')

case 5 %Z_out med sittig/bloomfield
        figure;

t=tiledlayout(2,1); % Her plottes outputtet af matricer som bodeplot
xlabel(t,'log spaced f values')

nexttile
semilogx(f,abs(Z_out)) % skriv enten F_DB, V_DB eller Z_DB
grid on
title('Magnitude plot of Z_{out}(f)')
ylabel('Magnitude of Z_{out}(dB)')

nexttile
semilogx(f,rad2deg(angle(Z_out))) % Samme her .
grid on
title('Angle plot of Z_{out}(f)')
ylabel('Angle of Z_{out}(^\circ)')
    otherwise
end

%bruh lad mig pull
