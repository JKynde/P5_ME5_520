%Kode til at sweepe forskellige parameter og lave et plot
clc
clear
j = sqrt(-1);
Parameters; % Først defineres alle parametre fra konstante Parameters.m
Tlmode = 1; % Tlmode til Matricer.m functionen. Der er om front layer matricen er ganget på. 1 for ja else ikke.
V_in = 300*exp(j*deg2rad(0)); % Indgangs spændingsvisor.
f=10*10^6;
% Lav en vektor logaritmisk fordelte indgange med n punkter.
N=1000;

c33D=zeros(N,1);
for n=(1:1:N)% (fra værdi : step size=1: til værdi)   
    c33D(n)=-N+2*n;
end
figure(1)
subplot(2,1,1)
plot(1:N,c33D)
xlabel('')
ylabel('')
grid