clc %hi hi
clear
Parameters;
f = 10^6; % Input frequency
omega = 2*pi*f; % Angular velocity
t=[0:1*10^-8*2:2*10^-6*2]; % initialize time vector t_start:t_step:t_end
z=[r_transducer*2:0.0001*0.5:r_transducer*2+2*lambda]; %initialize distance from transducer. The model is valid from 2*r_transducer
p_z=zeros(length(t),length(z)); %Pressure array.
Wavesum = zeros(length(t),length(z));% Initialize waves 
Wave1 = zeros(length(t),length(z));
Wave2 = zeros(length(t),length(z));
C=zeros(1,length(z)); %initialize diffraction coefficient.

F_zprvol=zeros(1,length(z)); % Vector for force pr. volume independent of time.

for i=1:length(t) % Loop over time and distance
    
for n=1:length(z)
[p_z(i,n),C(n),Wavesum(i,n),F_zprvol(n),Wave1(i,n),Wave2(i,n)]=Pressure(z(n),omega,t(i)); % Run pressure.m
end
Progress=(i/length(t))*100;
Progress
end

%% Post processesing and plotting


%plot(z,p_z(1,:),z,F_zprvol(:))
% 
% for i=1:length(t)
%    hold on
%    plot(z,Wavesum(i,:)); axis([2*r_transducer  z(length(z)) -2*10^6 2*10^6 ]); 
%    pause(0.1)
% end
    
    
    
    
for i=1:length(t)
   hold off
   plot(z,Wavesum(i,:),z,Wave1(i,:),z,Wave2(i,:)); axis([2*r_transducer  z(length(z)) -2*10^6 2*10^6 ]); 
   pause(0.1)
end
    









