clc %hi hi
clear
fprintf('Initializing \n')
Parameters;
%% Parametre og running options
plotte = 1; % plotting option 1 = plot mellem F_AC og P-felt til t =0 - 2 = plot af F_AC_ex og F_AC
generate_F_AC_ex_and_syms = 0; % Vælg om scriptet laver en function til kraften eller om den skal hoppe det over.
Tlmode=1; % Front lag eller ejjjj. Vi elsker voooores børn,
f=10^6; % Frekvens
V_in=300; % Spændingsfaser
omega=2*pi*f; % Vinkelhastighed
t_steps=50; % antal inddelinger i tid, skal bruges til num.int.
r_particle = 10^-4;
V_particle=(4/3)*pi*(r_particle)^3; %Volumen af partikel
z_stepsize=1.2500e-05/5; %Stepsize til z.
%% Initialisering
[F,v_t,~]=Matricer(f,V_in,Tlmode); %Henter F og v fra matricer
z=[r_transducer*2:z_stepsize:r_transducer*2+2*lambda]; %initaliserer afstandsvektoren fra 2 gange transducer radius, til 2 bølgelængder væk
t_stepsize=10^-6/(2*t_steps); %Størrelsen af steppet i tiden til simpsons rule i integral boy.
t=[0:t_stepsize:10^(-6)-t_stepsize]; % initialize time vector t_start:t_step:t_end t løber fra 0 til T (10^-6 s) minus en enkelt stepsize så det passer med at simpson får et lige antal punkter.
U_AC_V=zeros(1,length(z)); %Initialize primary radiation potential per volume vector 
F_AC_V=zeros(1,length(z)); %Initialize force per volume on a particle vector
F_AC=zeros(1,length(z)); %Initialize Force on a particle vector
WavesumresP=zeros(length(t),length(z)); %Laver et array med alle værdier for trykfeltet i både tid og afstand
Wavesumresv=zeros(length(t),length(z)); %Laver et array med alle værdier for hastighedsfeltet i både tid og afstand
WavesumresP_squared=zeros(length(t),length(z)); % Til p^2
Wavesumresv_squared=zeros(length(t),length(z)); % Til v^2

%% Long ass nested forloop
Progress=0;
fprintf('Filling wavesumres and squaring \n')
for n=1:length(z)
    for m=1:length(t)
    [~, ~, WavesumresP(m,n),~, ~, ~,Wavesumresv(m,n)]=Pressure(z(n),omega,t(m),F,v_t);  
    WavesumresP_squared(m,n)=(WavesumresP(m,n))^2;
    Wavesumresv_squared(m,n)=(Wavesumresv(m,n))^2;
    end
    Progress = (n/(length(z))*100);
    fprintf('Filling and squaring progress: %f \n',Progress)
end



% Wavesumresianden=(Wavesumres^2)
%% Integration og differentiation
P_avg=zeros(1,length(z));
v_avg=zeros(1,length(z));
fprintf('Integrating to find <p^2> and <v^2> \n')
for n=1:length(z)
   P_avg(n) = (1/(10^-6))*(integralboy(WavesumresP_squared(:,n),t_stepsize)); % trykket som funktion af tid og afstand i anden, integreret over tiden, divideret med perioden <p^2>
   v_avg(n) = (1/(10^-6))*(integralboy(Wavesumresv_squared(:,n),t_stepsize)); % hastigheden som funktion af tid og afstand i anden, integreret over tiden, divideret med perioden <v^2>
end
fprintf('Calculating U_AC_V \n')
for n=1:length(z)
   U_AC_V(n) = (f_2/(2*rho_oil*v_0Oil^2)*P_avg(n)-f_2*(3/4)*rho_oil*v_avg(n)); % Her udregnes gorkovs potential pr volumen af partikel
end

%% Plotting omkring F_AC
fprintf('Differentiating \n')
F_AC_V=-differentialboy(U_AC_V,z_stepsize);
F_AC = F_AC_V*V_particle;
if plotte == 1
hold on 
xlabel('Distance from transducer head')
yyaxis left
ylabel('Force on particle with given volume')
plot(z,F_AC)
yyaxis right 
ylabel('Pressurefield at time = 0')
plot(z,WavesumresP(1,:)) 
%plot (z,U_AC_V)
end

%% F_AC udtryksgenerator - Sinusbølge
if generate_F_AC_ex_and_syms == 1
    
% % F_AC function generator. F_AC antages at være en sinus funktion med
% given amplitude, bølgelængde og fasedrej, som findes ud fra den numeriske plot.

fprintf('Creating F_AC_function from numeric plot \n')

zeromap = zeros(1,3); % Laver vektor til n-værdier før fortegnskifte.
F_AC_max = max(F_AC); % Amplituden af F_AC findes.
% For loop to check for 3 sign flips (equal to three zeroes) and then
% lambda can be determined.
i = 1;

for n=1:length(F_AC) % Det her for loop går gennem F_AC vektoren og leder efter 
    % fortegnskift. Når den har fundet et fortegnskift på næste indgang,
    % gennem den tilsvarende n-værdi i zeromap og i tæller op. Når den har
    % fundet 3 nuller stopper den med at lede.

    if n ~= length(F_AC) && i <= 3
         if sign(F_AC(n)) ~= sign(F_AC(n+1))
            zeromap(i) = n;
            i=i+1;
         else
         end
    end
end

F_AC_lambda = z(zeromap(3))-z(zeromap(1)); % Bølgelængden udregnes, som forskellen mellem 
% det første og sidste nul, da man har 3 nuller in en enkelt bølgelængde.

F_AC_ex = F_AC_max*sin((2*pi/F_AC_lambda)*z); % Nu bygges der et udtryk for 
% F_AC (ex for expression), men den mangler fasedrejet. Bemærk at z
% vektoren indgør for at F_AC_ex udregnes for de samme z-værdier som F_AC.
% Derved passer n-værdierne og afstande og funktionsværdier kan
% sammenlignes.

zeromap2 = zeros(1,3);
i = 1;
for n=1:length(F_AC_ex) % Samme sang og dans gøres nu for F_AC_ex. Altså nullerne findes.

    if n ~= length(F_AC_ex) && i <= 3
         if sign(F_AC_ex(n)) ~= sign(F_AC_ex(n+1))
            zeromap2(i) = n;
            i=i+1;
         else
         end
    end
end
deltaZ=z(zeromap(1))-z(zeromap2(1)); % Nu kan vi finde afstanden mellem nullerne og 
% derved udrenge fasedrejet.
F_AC_ex = F_AC_max*sin((2*pi/F_AC_lambda)*z-2*pi/F_AC_lambda * deltaZ); % F_AC_ex opdateres
% med fasedrejet som 2pi/lambda * deltaZ.
syms x
F_AC_syms = F_AC_max*sin((2*pi/F_AC_lambda)*x-2*pi/F_AC_lambda * deltaZ); % Symbolic
% udtryk for F_AC er ovenstående. 
if plotte ==2
    
    
plot(z,F_AC_ex,z,F_AC); % Plot til at sammenligne F_AC_ex og F_AC
%vektorene
end

else
end


fprintf('Done \n')


