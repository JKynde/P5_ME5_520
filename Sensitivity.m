clc
clear
% Det her scribt tjekker for input parameter sensitivitet 
param=StructCreator();
OGparam=param;
fields = ["l_aa";"d_33";"l_a";"C_0"];
f = 10^6;
V_in = 150;
param.mode=3;
param.Tlmode = 1;
BaseCase = abs(Matricer2(f,V_in,param));
scalingfactor = 0.1;
upscale = 1+scalingfactor;%
downscale = 1-scalingfactor;%
format long

for n=1:length(fields)
fprintf("%s is being varied. Its value is originally %e\n",fields(n),param.(fields(n)))

param.(fields(n))=OGparam.(fields(n))*1.05;
fprintf("Upping %s with a factor %f. Its value is now %e\n",fields(n),upscale,param.(fields(n)))
Resultsup(n) = abs(Matricer2(f,V_in,param));
param.(fields(n))=OGparam.(fields(n))*0.95;
fprintf("Downing %s with a factor %f Its value is now %e \n",fields(n),downscale,param.(fields(n)))
Resultsdown(n)=abs(Matricer2(f,V_in,param));
param=OGparam;
end

for n=1:length(fields)
    bruh(n) = (Resultsup(n)-Resultsdown(n))/BaseCase ;
end
X = categorical(fields);
X = reordercats(X,fields);
bar(X,bruh)