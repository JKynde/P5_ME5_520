clc %hi hi
clear
Parameters;
f = 10^6; % Input frequency
V_in = 150; % Voltage phasor
omega = 2*pi*f; % Angular velocity
t=[0:1*10^-8*2:2*10^-6*2]; % initialize time vector t_start:t_step:t_end
z=[r_transducer*2:0.0001*0.5:r_transducer*2+2*lambda]; %initialize distance from transducer. The model is valid from 2*r_transducer
p_z=zeros(length(t),length(z)); %Pressure array.
Wavesum = zeros(length(t),length(z));% Initialize waves 
Wave1 = zeros(length(t),length(z));
Wave2 = zeros(length(t),length(z));
C=zeros(1,length(z)); %initialize diffraction coefficient.

F_zprvol=zeros(1,length(z)); % Vector for force pr. volume independent of time.

[F,v_t] = Matricer(omega/(2*pi),V_in,1); % Force and velocity from Sittig Model

for i=1:length(t) % Loop over time and distance
    
for n=1:length(z)
[p_z(i,n),C(n),Wavesum(i,n),F_zprvol(n),Wave1(i,n),Wave2(i,n),~]=Pressure(z(n),omega,t(i),F,v_t); % Run pressure.m
end
Progress=(i/length(t))*100;
Progress
end

%% Post processesing and plotting

% Plot over kraften af Lukes approximation sammen med trykket gange 2. Kraften samler sig i
% antinodesne pga. kontrastfaktoren phi (kan ses i trujillu er negativ.
%plot(z,F_zprvol(:))


% Plot over wavesum med enten hold on eller off. on = alle tidssteps kan
% ses samtidig. off = Der plottes hver tidsstep enkeltvis.
for i=1:length(t)
   hold on
   plot(z,Wavesum(i,:)); axis([2*r_transducer  z(length(z)) -2*10^6 2*10^6 ]); 
   pause(0.1)
end

    
%     % Plot over Wave1, Wave 2 og Wavesum som de flytter sig gennem space.
% for i=1:length(t)
%    hold off
%    plot(z,Wavesum(i,:),z,real(Wave1(i,:)),z,real(Wave2(i,:))); axis([2*r_transducer  z(length(z)) -2*10^6 2*10^6 ]); 
%    pause(0.1)
% end
%     









