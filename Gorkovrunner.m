clc %hi hi
clear
Parameters;
Tlmode=1; % Front lag eller ejjjj. Vi elsker voooores børn,
f=10^6; % Frekvens
V_in=300; % Spændingsfaser
omega=2*pi*f; % Vinkelhastighed
t_steps=50; % antal inddelinger i tid, skal bruges til num.int.
r_particle = 10^-4;
V_particle=(4/3)*pi*(r_particle)^3; %Volumen af partikel
z_stepsize=0.0001*0.5/4;

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

for n=1:length(z)
    for m=1:length(t)
    [~, ~, WavesumresP(m,n),~, ~, ~,Wavesumresv(m,n)]=Pressure(z(n),omega,t(m),F,v_t);  
    WavesumresP_squared(m,n)=(WavesumresP(m,n))^2;
    Wavesumresv_squared(m,n)=(Wavesumresv(m,n))^2;
    end
end



% Wavesumresianden=(Wavesumres^2)

P_avg=zeros(1,length(z));
v_avg=zeros(1,length(z));

for n=1:length(z)
   P_avg(n) = (1/(10^-6))*(integralboy(WavesumresP_squared(:,n),t_stepsize)); % trykket som funktion af tid og afstand i anden, integreret over tiden, divideret med perioden <p^2>
   v_avg(n) = (1/(10^-6))*(integralboy(Wavesumresv_squared(:,n),t_stepsize)); % hastigheden som funktion af tid og afstand i anden, integreret over tiden, divideret med perioden <v^2>
end

for n=1:length(z)
   U_AC_V(n) = (f_2/(2*rho_oil*v_0Oil^2)*P_avg(n)-f_2*(3/4)*rho_oil*v_avg(n)); % Her udregnes gorkovs potential pr volumen af partikel
end

F_AC_V=-differentialboy(U_AC_V,z_stepsize);
F_AC = F_AC_V*V_particle;
hold on 
xlabel('Distance from transducer head')
yyaxis left
ylabel('Force on particle with given volume')
plot(z,F_AC)
yyaxis right 
ylabel('Pressurefield at time = 0')
plot(z,WavesumresP(1,:))
