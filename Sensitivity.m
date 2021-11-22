clc
clear
% Det her scribt tjekker for input parameter sensitivitet 
param=StructCreator();
fields = fieldnames(param);
f = 10^6;
V_in = 150;
param.mode=2;
param.Tlmode = 1;
BaseCase = Matricer2(f,V_in,param);

boi="C_0"; % boi er den parameter, som man godt kunne tænke sig at undersøge sensitiviteten af.
param.(boi)
