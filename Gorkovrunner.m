clc %hi hi
clear
Parameters;
Tlmode=1; % Front lag eller ejjjj. Vi elsker voooores børn,
f=10^6; % Frekvens
V_in=300; % Spændingsfaser
omega=2*pi*f; % Vinkelhastighed
t_steps=50; % antal inddelinger i tid, skal bruges til num.int.
r_particle = 10^-4;
V_particle=(4/3)*pi*(r_particle)^3; %Volumen af partikel, bruges i gorkov
z_stepsize=0.0001*0.5;


[F,v_t,~]=Matricer(f,V_in,Tlmode); %Henter F og v fra matricer
z=[r_transducer*2:z_stepsize:r_transducer*2+2*lambda]; %initaliserer afstandsvektoren fra 2 gange transducer radius, til 2 bølgelængder væk
t_stepsize=10^-6/(2*t_steps); %Størrelsen af steppet i tiden til simpsons rule i integral boy.
t=[0:t_stepsize:10^(-6)-t_stepsize]; % initialize time vector t_start:t_step:t_end
U_AC_V=zeros(length(z)); %Initialize primary radiation potential per volume vector 
F_AC_V=zeros(length(z)); %Initialize force per volume on a particle vector
F_AC=zeros(length(z)); %Initialize Force on a particle vector
WavesumresP=zeros(length(t),length(z)); %Laver et array med alle værdier for trykfeltet i både tid og afstand
Wavesumresv=zeros(length(t),length(z)); %Laver et array med alle værdier for hastighedsfeltet i både tid og afstand
for n=1:length(z)
    for m=1:length(t)
    [~, ~, WavesumresP(m,n),~, ~, ~,Wavesumresv(m,n)]=Pressure(z(n),omega,t(m),F,v_t);   
    end
end
% Wavesumresianden=(Wavesumres^2)
P_avg=zeros(1,length(z));
v_avg=zeros(1,length(z));
for n=1:length(z)
   P_avg(n) = (1/(10^-6)*integralboy(WavesumresP(:,n),t_stepsize))^2; % trykket som funktion af tid og afstand i anden, integreret over tiden, divideret med perioden <p^2>
   v_avg(n) = (1/(10^-6)*integralboy(Wavesumresv(:,n),t_stepsize))^2; % hastigheden som funktion af tid og afstand i anden, integreret over tiden, divideret med perioden <v^2>
end


U_AC_V =zeros(1,length(z)) ; % force potential pr. volume 

for n=1:length(z)
   U_AC_V(n) = (f_2/(2*rho_oil*v_0Oil^2)*P_avg(n)-f_2*(3/4)*rho_oil*v_avg(n)); % Her udregnes gorkovs potential pr volumen af partikel
end

F_AC_V=differentialboy(U_AC_V,z_stepsize);
F_AC = F_AC_V*V_particle;
hold on 
xlabel('Distance from transducer head')
yyaxis left
ylabel('Force on particle with given volume')
plot(z,F_AC)
yyaxis right 
ylabel('Pressurefield at time = 0')
plot(z,WavesumresP(1,:))
