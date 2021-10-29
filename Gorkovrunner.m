clc %hi hi
clear
Parameters;
Tlmode=1;
f=10^6;
V_in=300;
omega=2*pi*f;
t_steps=50; % antal inddelinger i tid, skal bruges til num.int.
[F,v_t,~]=Matricer(f,V_in,Tlmode); %Henter
z=[r_transducer*2:0.0001*0.5:r_transducer*2+2*lambda]; %initaliserer afstandsvektoren fra 2 gange transducer radius, til 2 bølgelængder væk
t_stepsize=10^-6/(2*t_steps); %Størrelsen af steppet i tiden
t=[0:t_stepsize:10^(-6)-t_stepsize]; % initialize time vector t_start:t_step:t_end
U_AC_V=zeros(length(z)); %Initialize primary radiation potential per volume vector 
F_AC_V=zeros(length(z)); %Initialize force per volume on a particle vector
F_AC=zeros(length(z)); %Initialize Force on a particle vector
Wavesumres=zeros(length(t),length(z)); %Laver et array med alle værdier for trykfeltet i både tid og afstand
for n=1:length(z)
    for m=1:length(t)
    [~, ~, Wavesumres(m,n),~, ~, ~]=Pressure(z(n),omega,t(m),F,v_t);   
    end
end
% Wavesumresianden=(Wavesumres^2)
P_avg=zeros(1,length(z));
% for n=1:length(z)
%    P_avg=   % trykket som funktion af tid og afstand i anden, integreret over tiden, divideret med perioden  
% end
   
