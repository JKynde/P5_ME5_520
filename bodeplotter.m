clc
clear all
param=StructCreator();
OGparam=param;
f = logspace(5,6.5,5000);
V_in = 150;
param.mode = 3;
param.Tlmode = 1;
j = sqrt(-1);
A = zeros(1,length(f));
B = zeros(1,length(f));
C = zeros(1,length(f));
for n = 1:length(f)
    [~,~,~,A(n)] = Matricer2(f(n),V_in,param);
end
param.rho_P=param.rho_p*0.5;
param = StructRecalculator(param)
for n = 1:length(f)
    [~,~,~,B(n)] = Matricer2(f(n),V_in,param);
end
param = OGparam;
param.c33D=param.c33D*10;
param = StructRecalculator(param); 
for n = 1:length(f)
    [~,~,~,C(n)] = Matricer2(f(n),V_in,param);
end
A = 20*log10(abs(A));
B = 20*log10(abs(B));
C = 20*log10(abs(C));
figure;
semilogx(f, A, f, B, f, C)
grid on
legend({'original', 'scaled up', 'scaled down'}, 'FontSize', 14, 'Location', 'northwest');
title('Amplitude of frequency response', 'FontSize', 20);
xlabel('Logarithmically spaced frequency values', 'FontSize', 20);
ylabel('dB', 'FontSize', 20);